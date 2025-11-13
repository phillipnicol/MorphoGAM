#' @export
#'
#' @title CurveFinder: Automatically estimate curve in two-dimensional ST data
#'
#' @description This function estimates a (possibly closed) parametric curve
#' passing through spatial transcriptomics data and uses this curve to define
#' morphologically relevant coordinates `t` and `r`.
#'
#' @param xy A numeric matrix or data frame with two columns representing the x and y coordinates of the data points.
#' @param knn An integer specifying the number of nearest neighbors used to construct the KNN graph. Default is 5.
#' @param prune.outlier A numeric threshold for pruning outliers based on the distance to the k+1 nearest neighbor. Outliers are removed if their distance exceeds `prune.outlier * median(nnk)`. Defaults to NULL (no pruning).
#' @param loop A logical value indicating whether the curve should be treated as a loop (closed curve). Default is FALSE.
#'
#' @return A list containing:
#' \item{xyt}{A data frame with x and y coordinates, fitted curve parameters (`t` and `r`), and fitted values for x and y (`f1` abd `f2`).}
#' \item{curve.plot}{A ggplot object showing the original data with the fitted curve overlaid.}
#' \item{coordinate.plot}{A ggplot object displaying the data points color-coded by their fitted `t` values.}
#' \item{residuals.plot}{A ggplot object displaying the data points color-coded by their fitted `r` values.}
#'
#' @details
#' To-do
#'
#'
#' @references
#' Kraemer G, Reichstein M, Mahecha MD (2018). “dimRed and coRanking—Unifying Dimensionality Reduction in R.” _The R Journal_, *10*(1), 342-358.
#'
#'
#' @import ggplot2
#' @import mgcv
#' @importFrom igraph components subgraph V
#' @importFrom RSpectra eigs_sym
#' @importFrom gratia derivatives
#' @importFrom gtools permutations
#' @importFrom RANN nn2
CurveFinder <- function(xy,
                       knn=5,
                       prune.outlier=NULL,
                       loop=FALSE) {
  xy.dist <- as.matrix(dist(xy))

  if(!is.null(prune.outlier)) {
    nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
    outlier <- which(nnk > (prune.outlier)*median(nnk))
    if(length(outlier) > 0) {
      xy.outlier <- xy[outlier,]
      xy <- xy[-outlier,]
    }
  }

  knng <- dimRed:::makeKNNgraph(x = xy,
                                k = knn,
                                eps = 0)

  if (!is.logical(loop) && loop != "auto") {
    stop("'loop' must be one of TRUE, FALSE, or \"auto\"")
  } else if (loop == "auto") {
    loop <- is_loop(xy, knng)
    message(ifelse(loop, "Detected loop.", "Detected path."))
  }

  comp <- igraph::components(knng)

  if(comp$no > 5 || (loop && comp$no > 1)) {
    stop("Graph has too many disconnected components. Increase knn")
  }

  endpoints <- matrix(0, nrow=2*comp$no, ncol=2)
  t.list <- list()
  for(c in 1:comp$no) {
    xy.sub <- xy[comp$membership == c,]

    knng.sub <- igraph::subgraph(knng,
                                 igraph::V(knng)[comp$membership == c])

    # See ISOMAP function in
    ##  Kraemer G, Reichstein M, Mahecha MD (2018). “dimRed
    #and coRanking-Unifying Dimensionality Reduction in
    #R.” _The R Journal_, *10*(1), 342-358. coRanking
    #version 0.2.6,
    #<https://journal.r-project.org/archive/2018/RJ-2018-039/index.html>.

    dG <- igraph::distances(knng.sub, algorithm = "dijkstra")

    dG <- dG ^ 2
    dG <- .Call(stats:::C_DoubleCentre, dG)
    dG <- - dG/2

    if(loop) {
      e <- RSpectra::eigs_sym(dG, 2, which = "LA",
                              opts = list(retvec = TRUE))
      t <- (pi+atan2(e$vectors[,1], e$vectors[,2]))/(2*pi)
    } else{
      e <- RSpectra::eigs_sym(dG, 1, which = "LA",
                              opts = list(retvec = TRUE))

      t.list[[c]] <- as.vector(e$vectors)

      t.list[[c]] <- (t.list[[c]] - min(t.list[[c]]))/(max(t.list[[c]]) - min(t.list[[c]]))
      t.list[[c]] <- t.list[[c]]*(sum(comp$membership == c)/length(comp$membership))

      endpoints[2*c-1,] <- xy.sub[which.min(t.list[[c]]),]
      endpoints[2*c,] <- xy.sub[which.max(t.list[[c]]),]
    }
  }

  if(comp$no > 1) {
    paths <- segment_TSP(endpoints)
    best.path <- which.min(paths$cost)
    t.shift <- 0
    t <- rep(0, nrow(xy))
    for(c in 1:comp$no) {
      current.comp <- paths[best.path,c]

      if(current.comp > 0) {
        t[comp$membership == current.comp] <- t.list[[current.comp]] + t.shift
      } else{
        t.max <- max(t.list[[-1*current.comp]])
        t[comp$membership == -1*current.comp] <- t.max - t.list[[-1*current.comp]] + t.shift
      }
      t.shift <- max(t)
    }
  } else if(!loop){
    t <- as.vector(t.list[[1]])
  }

  knots <- round(min(100, 0.1*nrow(xy)))
  if(loop) {
    fitx <- mgcv::gam(xy[,1]~s(t,bs="cc",k=knots))
    fity <- mgcv::gam(xy[,2]~s(t,bs="cc",k=knots))
  } else{
    fitx <- mgcv::gam(xy[,1]~s(t,bs="cr",k=knots))
    fity <- mgcv::gam(xy[,2]~s(t,bs="cr",k=knots))
  }

  r <- orthogonal_path(fitx,fity,t)

  my.t <- seq(0,1,by=10^{-4})
  predx <- predict(fitx,newdata=list(t=my.t))
  predy <- predict(fity,newdata=list(t=my.t))
  ft <- data.frame(x=predx,y=predy)

  df.new <- data.frame(x=xy[,1],
                       y=xy[,2])

  p <- ggplot(data=df.new,aes(x=x,y=y)) + geom_point(col="grey",
                                                     alpha=0.5)

  df.line <- data.frame(x=ft[,1], y=ft[,2], color=my.t)
  p <- p + geom_path(data=df.line,aes(x=x,y=y,color=color),
                     linewidth=1)
  p <- p + scale_color_gradient(low="navyblue",
                                high="firebrick1") +
    labs(color="t")
  p <- p + theme_bw()

  xyt <- data.frame(x=xy[,1],y=xy[,2],t=t, r=r,
                    f1 = fitted(fitx),
                    f2=fitted(fity))


  p2 <- data.frame(x=xy[,1],y=xy[,2],color=t) |>
    ggplot(aes(x=x,y=y,color=color)) + geom_point() +
    scale_color_gradientn(values=c(0,0.5,1),
                          colors=c("navyblue","grey90", "firebrick1"))+
    theme_bw() +
    ggtitle("First Coordinate") +
    labs(color="t")

  p3 <- data.frame(x=xy[,1],y=xy[,2],color=r) |>
    ggplot(aes(x=x,y=y,color=color)) + geom_point() +
    scale_color_gradientn(values=c(0,0.5,1),
                          colors=c("navyblue","grey90", "firebrick1"))+
    theme_bw() +
    ggtitle("Second Coordinate") +
    labs(color="r")

  out <- list()
  out$xyt <- xyt
  out$curve.plot <- p
  out$coordinate.plot <- p2
  out$residuals.plot <- p3
  return(out)
}

orthogonal_path <- function(fitx,fity, t) {
  f2x <- gratia::derivatives(fitx,order=1,data=data.frame(t=t))
  f2y <- gratia::derivatives(fity,order=1,data=data.frame(t=t))
  t2 <- t

  for(i in 1:length(t)) {
    e <- c(fitx$residuals[i], fity$residuals[i])
    Rf1 <- c(-f2y$.derivative[i], f2x$.derivative[i])
    sign <- ifelse(sum(e*Rf1) > 0, 1, -1)
    #t2[i] <- sum(e[1]*Rf1[1] + e[2]*Rf1[2])
    t2 <- sign*sqrt(sum(e^2))
  }

  t2 <- (t2 - min(t2))/(max(t2) - min(t2))
  return(t2)
}


segment_TSP <- function(endpoints) {
  ncomp <- nrow(endpoints)/2

  orientation <- expand.grid(replicate(ncomp, c(-1,1), simplify = FALSE))

  paths <- matrix(0, nrow=0,ncol=ncomp)
  for(i in 1:nrow(orientation)) {
    vec <- 1:ncomp
    vec[orientation[i,] == -1] <- -1*vec[orientation[i,] == -1]
    paths <- rbind(paths, gtools::permutations(n=ncomp, r=ncomp, v=vec))
  }

  my.dist <- as.matrix(dist(endpoints))
  costs <- apply(paths, 1, function(x) {
    cost <- 0
    for(j in 1:(ncomp - 1)) {
      if(x[j] < 0) {
        if(x[j+1] < 0) {
          cost <- cost + my.dist[-2*x[j]-1, -2*x[j+1]]
        } else{
          cost <- cost + my.dist[-2*x[j]-1, 2*x[j+1]-1]
        }
      } else{
        if(x[j+1] < 0) {
          cost <- cost + my.dist[2*x[j], -2*x[j+1]]
        } else{
          cost <- cost + my.dist[2*x[j], 2*x[j+1]-1]
        }
      }
    }
    cost
  })

  colnames(paths) <- 1:ncomp
  paths <- as.data.frame(paths)
  paths$cost <- as.vector(costs)

  return(paths)
}


### Taken from package dimRed due to it being delisted from CRAN
#   Kraemer G, Reichstein M, Mahecha MD (2018). “dimRed and
#coRanking-Unifying Dimensionality Reduction in R.” _The R
#Journal_, *10*(1), 342-358. coRanking version 0.2.6,
#<https://journal.r-project.org/archive/2018/RJ-2018-039/index.html>.
dimRed_makeKNNgraph <- function (x, k, eps = 0, diag = FALSE)
{
  INF_VAL <- 1.340781e+15
  NA_IDX <- 0
  BDKD_LIM <- 1e+06
  M <- nrow(x)
  treetype <- "kd"
  searchtype <- if (eps == 0)
    "standard"
  else "priority"
  nn2res <- RANN::nn2(data = x, query = x, k = k + 1, treetype = treetype,
                      searchtype = searchtype, eps = eps)
  g <- igraph::make_empty_graph(M, directed = TRUE)
  g[from = if (diag)
    rep(seq_len(M), times = k + 1)
    else rep(seq_len(M), times = k), to = if (diag)
      as.vector(nn2res$nn.idx)
    else as.vector(nn2res$nn.idx[, -1]), attr = "weight"] <- if (diag)
      as.vector(nn2res$nn.dists)
  else as.vector(nn2res$nn.dists[, -1])
  return(igraph::as.undirected(g, mode = "collapse", edge.attr.comb = "first"))
}







