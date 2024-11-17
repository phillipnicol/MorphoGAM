

#' @export
#'
#' #' @title scDist: Identify perturbed cell types in
#' single-cell RNA-seq data
#'
#' @description Estimate the distance between
#' condition means in gene expression space.
#'
#' @param normalized_counts A matrix containing
#' normalized data with genes on rows and cells
#' on columns
#' @param meta.data A data frame containing meta data for each cell.
#' @param fixed.effects The columns in meta.data corresponding to the fixed effects. In
#' a typical case, this would be the condition of interest.
#' @param random.effects The columns in meta.data corresponding to the random effects.
#' In a typical use case this would be the column containing patient ID.
#' @param clusters The column containing the cell-type annotation.
#' @param d The number of PCs to use.
#' @param truncate Whether or not to round negative distances to 0.
#' @param min.count.per.cell The minimum number of cells per cluster to perform the estimation.
#' @param weights An optional vector of length equal to the number of genes specifying the weight
#' to place on each gene in the distance estimate.
#'
#' @return A list with components
#' \itemize{
#' \item \code{results} - A data frame containing the cell
#' type, estimated distance, and other statistics such as p-value.
#' \item \code{vals} For each cell type a list of more detailed
#' information (such as raw data) and coefficients for each PC are
#' included.
#' }
#'
#'
#' @title MorphoGAM
#'
#' @importFrom ashr ash get_post_sample
MorphoGAM <- function(Y,
                      curve.fit,
                      design,
                      shrinkage=TRUE,
                      min.count.per.gene=10,
                      return.fx = TRUE,
                      offset=NULL) {

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
    fit <- mgcv::gam(formula=design,
                     family=mgcv::nb(),
                     offset=l.o,
                     data=data,
                     H=H)

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
    if(length(t.cols) > 0) {
      fx.t <- basis.functions[,t.cols] %*% beta.shrink[t.cols]
      peak.t <- max(abs(fx.t))
      range.t <- max(exp(beta_g0 + fx.t)) - min(exp(beta_g0 + fx.t))
      range.t <- median.depth*range.t
      pv.t <- summary(fit)$s.table["s(t)",4]
      results[i,c("peak.t","range.t","pv.t")] <- c(peak.t,range.t,pv.t)
      if(return.fx) {
        fxs.t[i,] <- fx.t
      }
    }

    r.cols <- grep("s\\(r\\)", colnames(basis.functions))
    if(length(r.cols) > 0) {
      fx.r <- basis.functions[,r.cols] %*% beta.shrink[r.cols]
      peak.r <- max(abs(fx.r))
      range.r <- max(exp(beta_g0 + fx.r)) - min(exp(beta_g0 + fx.r))
      range.r <- median.depth*range.r
      pv.r <- summary(fit)$s.table["s(r)",4]
      results[i,c("peak.r","range.r","pv.r")] <- c(peak.r,range.r,pv.r)
      if(return.fx) {
        fxs.r[i,] <- fx.r
      }
    }
  }

  out <- list()
  out$results <- as.data.frame(results)
  rownames(out$results) <- rownames(Y)

  if(return.fx) {
    rownames(fxs.t) <- rownames(Y)
    rownames(fxs.r) <- rownames(Y)
    out$fxs.t <- fxs.t
    out$fxs.r <- fxs.r
  }

  return(out)
}

