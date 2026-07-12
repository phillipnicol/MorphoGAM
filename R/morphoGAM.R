
.add_default_cr_smooth_basis <- function(expr) {
  if(!is.call(expr)) {
    return(expr)
  }

  call.name <- expr[[1]]
  is.smooth <- identical(call.name, as.name("s")) ||
    (is.call(call.name) &&
       identical(call.name[[1]], as.name("::")) &&
       identical(call.name[[3]], as.name("s")))

  if(is.smooth) {
    if(!("bs" %in% names(as.list(expr)))) {
      expr$bs <- "cr"
    }
    return(expr)
  }

  for(i in seq_along(expr)[-1]) {
    expr[[i]] <- .add_default_cr_smooth_basis(expr[[i]])
  }

  expr
}

.default_cr_smooth_basis <- function(design) {
  design <- stats::as.formula(design)

  for(i in seq_along(design)[-1]) {
    design[[i]] <- .add_default_cr_smooth_basis(design[[i]])
  }

  design
}

#' @export
#'
#' @title MorphoGAM: Apply a GAM to morphologically relevant coordinates
#'
#' @description
#' Fit gene-wise negative-binomial GAMs to identify expression variation with
#' respect to morphologically relevant coordinates estimated by `CurveFinder()`
#' or `CurveFinderInteractive()`.
#'
#' @param Y A numeric count matrix with genes in rows and samples or cells in
#'   columns.
#' @param curve.fit An object produced by `CurveFinder()` or
#'   `CurveFinderInteractive()`, containing fitted `t` and `r` coordinates for
#'   each sample.
#' @param design A formula specifying the GAM design, typically including smooth
#'   terms for `t` and `r` such as `y ~ s(t) + s(r)`. Smooth terms that do not
#'   specify `bs` use `bs = "cr"` by default.
#' @param shrinkage Logical; if `TRUE`, shrink smooth term coefficients using
#'   `ashr`. Default is `FALSE`.
#' @param min.count.per.gene Minimum total count required for a gene to be fit.
#'   Genes below this threshold are returned with null/default statistics.
#' @param return.fx Logical; if `TRUE`, return fitted smooth terms and FPCA
#'   summaries for the `t` and `r` effects.
#' @param offset Optional numeric vector of model offsets. If `NULL`, offsets
#'   are computed as the logarithm of column sums of `Y`.
#' @param knots.t Optional numeric vector specifying knot locations for the
#'   smooth term in `t`.
#' @param knots.r Optional numeric vector specifying knot locations for the
#'   smooth term in `r`.
#'
#' @return A list containing:
#' \item{results}{A data frame with one row per gene and columns `peak.t`,
#' `range.t`, `pv.t`, `peak.r`, `range.r`, `pv.r`, and `intercept`.}
#' \item{design}{The model formula used after adding default cubic regression
#' spline bases to unspecified smooth terms.}
#' If `return.fx = TRUE`, the list also includes:
#' \item{fxs.t}{A matrix of fitted smooth terms for `t` (genes by samples).}
#' \item{fxs.r}{A matrix of fitted smooth terms for `r` (genes by samples).}
#' \item{fpca.t}{An `irlba` decomposition of `fxs.t`, when available.}
#' \item{fpca.r}{An `irlba` decomposition of `fxs.r`, when available.}
#'
#' @importFrom ashr ash get_post_sample
#' @importFrom irlba irlba
#' @importFrom stats median predict vcov
#' @importFrom utils setTxtProgressBar txtProgressBar
MorphoGAM <- function(Y,
                      curve.fit,
                      design,
                      shrinkage=FALSE,
                      min.count.per.gene=10,
                      return.fx = TRUE,
                      offset=NULL,
                      knots.t = NULL,
                      knots.r = NULL) {

  design <- .default_cr_smooth_basis(design)

  if(is.null(offset)) {
    l.o <- log(colSums(Y))
  } else{
    #No offset
    l.o <- offset
  }
  median.depth <- median(colSums(Y))
  n <- ncol(Y); p <- nrow(Y)

  data <- data.frame(y = Y[1,],
                     t = curve.fit$xyt$t,
                     r = curve.fit$xyt$r,
                     l.o=l.o)

  fit.setup <- mgcv::gam(formula=design,family=mgcv::nb(),
                         data=data, fit=FALSE)

  H <- diag(ncol(fit.setup$X)); H[1,1] <- 0 #Don't penalize intercept

  #Setup results data frame
  results <- matrix(0, nrow=p, ncol=7)
  colnames(results) <- c("peak.t", "range.t", "pv.t",
                         "peak.r", "range.r", "pv.r", "intercept")

  if(return.fx) {
    fxs.t <- matrix(0,nrow=p,ncol=n)
    fxs.r <- matrix(0,nrow=p,ncol=n)
  }


  bar <- txtProgressBar(min=0,max=nrow(Y)-1,initial = 0)
  for(i in 1:nrow(Y)) {
    setTxtProgressBar(bar,i)
    if(rowSums(Y)[i] < min.count.per.gene) {
      results[i,] <- c(0,0,1,0,0,1,0)
      next
    }

    data$y <- Y[i,]

    smooth.knots <- list()
    if(!is.null(knots.t)) {
      smooth.knots$t <- knots.t
    }
    if(!is.null(knots.r)) {
      smooth.knots$r <- knots.r
    }

    if(length(smooth.knots) == 0) {
      fit <- mgcv::gam(formula=design,
                       family=mgcv::nb(),
                       offset=l.o,
                       data=data,
                       H=H)
    } else{
      fit <- mgcv::gam(formula=design,
                       family=mgcv::nb(),
                       offset=l.o,
                       data=data,
                       H=H,
                       knots=smooth.knots)
    }


    beta_g0 <- fit$coefficients[1]
    results[i,"intercept"] <- beta_g0

    se_beta <- diag(vcov(fit))[-1] #remove intercept

    #Shrink
    if(shrinkage) {
      my.ash <- ashr::ash(fit$coefficients[-1], se_beta)
      beta.shrink <- apply(ashr::get_post_sample(my.ash,1000),
                           2,
                           median)
    } else {
      beta.shrink <- fit$coefficients[-1]
    }

    basis.functions <- mgcv::predict.gam(fit, type="lpmatrix")[,-1]



    t.cols <- grep("s\\(t\\)", colnames(basis.functions))

    pred <- predict(fit, newdata = data.frame(t=data$t,r=data$r),
                    type="terms", se.fit=TRUE)
    if(length(t.cols) > 0) {
      fx.t <- basis.functions[,t.cols] %*% beta.shrink[t.cols]

      fx.t <- pred$fit[,1]/(1 + pred$se.fit[,1]^2)

      if(return.fx) {
        fxs.t[i,] <- fx.t
      }

      fx.t <- ifelse(pred$fit[,1] > 1.96*pred$se.fit[,1],
                     pred$fit[,1] - 1.96*pred$se.fit[,1],
                     ifelse(pred$fit[,1] < -1.96*pred$se.fit[,1],
                            pred$fit[,1] + 1.96*pred$se.fit[,1],
                            0))
      peak.t <- max(abs(fx.t))
      range.t <- max(exp(beta_g0 + fx.t)) - min(exp(beta_g0 + fx.t))
      range.t <- median.depth*range.t
      pv.t <- summary(fit)$s.table["s(t)",4]
      results[i,c("peak.t","range.t","pv.t")] <- c(peak.t,range.t,pv.t)
    }

    r.cols <- grep("s\\(r\\)", colnames(basis.functions))
    if(length(r.cols) > 0) {
      fx.r <- basis.functions[,r.cols] %*% beta.shrink[r.cols]
      fx.r <- pred$fit[,2]/(1 + pred$se.fit[,2]^2)

      if(return.fx) {
        fxs.r[i,] <- fx.r
      }

      fx.r <- ifelse(pred$fit[,2] > 1.96*pred$se.fit[,2],
                        pred$fit[,2] - 1.96*pred$se.fit[,2],
                        ifelse(pred$fit[,2] < -1.96*pred$se.fit[,2],
                               pred$fit[,2] + 1.96*pred$se.fit[,2],
                               0))

      peak.r <- max(abs(fx.r))
      range.r <- max(exp(beta_g0 + fx.r)) - min(exp(beta_g0 + fx.r))
      range.r <- median.depth*range.r
      pv.r <- summary(fit)$s.table["s(r)",4]
      results[i,c("peak.r","range.r","pv.r")] <- c(peak.r,range.r,pv.r)
    }
  }

  out <- list()
  out$results <- as.data.frame(results)
  rownames(out$results) <- rownames(Y)

  if(return.fx) {
    rownames(fxs.t) <- rownames(Y)
    rownames(fxs.r) <- rownames(Y)
    try({
      out$fpca.t <- irlba::irlba(fxs.t)
      out$fpca.r <- irlba::irlba(fxs.r)
    })
    out$fxs.t <- fxs.t
    out$fxs.r <- fxs.r
  }

  #Also return design
  out$design <- design
  return(out)
}
