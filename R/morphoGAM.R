

#' @export
#'
#' @title MorphoGAM: Apply a GAM to the morphologically relevant coordinates
#'
#' @description The morphologically relevant coordinates can be estimated using
#' `CurveFinder()` or `CurveFinderInteractive()`. This function applies a flexible
#' count-based model to identify genes with spatially variable expression with
#' respect to the morphologically relevant coordinates.
#'
#' @description Estimate the distance between
#' condition means in gene expression space.
#'
#' @param Y A numeric matrix where rows represent genes and columns represent samples (e.g., cells).
#' @param curve.fit An object produced by `CurveFinder`, containing the fitted curve parameters (`t` and `r`) for each sample.
#' @param design A formula specifying the GAM design, typically including smooth terms for `t` and `r` (e.g., `y ~ s(t) + s(r)`).
#' @param shrinkage A logical value indicating whether to apply shrinkage to smooth term coefficients using `ashr`. Default is `TRUE`.
#' @param min.count.per.gene An integer specifying the minimum total count required for a gene to be included in the analysis. Default is 10.
#' @param return.fx A logical value indicating whether to return the fitted smooth terms for `t` and `r` for each gene. Default is `TRUE`.
#' @param offset A numeric vector providing offset values for the GAM model. If `NULL` (default), offsets are computed as the logarithm of the column sums of `Y`.
#'
#' @return A list containing:
#' \item{results}{A data frame with rows corresponding to genes and the following columns:
#' \describe{
#'   \item{peak.t}{Maximum absolute smooth term for `t`.}
#'   \item{range.t}{Range of predicted expression values along `t`.}
#'   \item{pv.t}{p-value for the smooth term `t`.}
#'   \item{peak.r}{Maximum absolute smooth term for `r`.}
#'   \item{range.r}{Range of predicted expression values along `r`.}
#'   \item{pv.r}{p-value for the smooth term `r`.}
#'   \item{intercept}{Intercept of the fitted model.}
#' }}
#' If `return.fx = TRUE`, the list also includes:
#' \item{fxs.t}{A matrix of fitted smooth terms for `t` (genes x samples).}
#' \item{fxs.r}{A matrix of fitted smooth terms for `r` (genes x samples).}
#'
#'
#' @importFrom ashr ash get_post_sample
#' @importFrom irlba irlba
MorphoGAM <- function(Y,
                      curve.fit,
                      design,
                      shrinkage=FALSE,
                      min.count.per.gene=10,
                      return.fx = TRUE,
                      offset=NULL,
                      knots.t = NULL,
                      knots.r = NULL) {

  if(is.null(offset)) {
    l.o <- log(colSums(Y))
  } else{
    l.o <- rep(0,ncol(Y))
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

    if(is.null(knots.t)) {
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
                       knots=list(t = knots.t))
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

      # Alternative shrinkage
      #fx.t <- ifelse(pred$fit[,1] > 1.96*pred$se.fit[,1],
      #               pred$fit[,1] - 1.96*pred$se.fit[,1],
      #               ifelse(pred$fit[,1] < -1.96*pred$se.fit[,1],
      #                      pred$fit[,1] + 1.96*pred$se.fit[,1],
      #                      0))
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

      #fx.r <- ifelse(pred$fit[,2] > 1.96*pred$se.fit[,2],
      #                  pred$fit[,2] - 1.96*pred$se.fit[,2],
      #                  ifelse(pred$fit[,2] < -1.96*pred$se.fit[,2],
      #                         pred$fit[,2] + 1.96*pred$se.fit[,2],
      #                         0))

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

  return(out)
}

