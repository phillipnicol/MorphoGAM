

CurveFinder <- function(xy,
                       knn=5,
                       prune.outlier=NULL) {
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

  comp <- components(knng)

  if(comp$no > 5) {
    stop("Graph has too many disconnected components. \n
         Increase knn")
  }

  endpoints <- matrix(0, nrow=2*comp$no, ncol=2)
  t.list <- list()
  for(c in 1:comp$no) {
    xy.sub <- xy[comp$membership == c,]

    knng.sub <- subgraph(knng, V(knng)[comp$membership == c])

    geodist <- igraph::distances(knng.sub, algorithm = "dijkstra")

    k <- geodist ^ 2
    k <- .Call(stats:::C_DoubleCentre, k)
    k <- - k / 2

    e <- RSpectra::eigs_sym(k, 1, which = "LA",
                            opts = list(retvec = TRUE))

    t.list[[c]] <- as.vector(e$vectors)

    t.list[[c]] <- (t.list[[c]] - min(t.list[[c]]))/(max(t.list[[c]]) - min(t.list[[c]]))
    t.list[[c]] <- t.list[[c]]*(sum(comp$membership == c)/length(comp$membership))

    endpoints[2*c-1,] <- xy.sub[which.min(t.list[[c]]),]
    endpoints[2*c,] <- xy.sub[which.max(t.list[[c]]),]
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
  } else{
    t <- as.vector(t.list[[1]])
  }

  knots <- round(min(100, 0.1*nrow(xy)))
  fitx <- gam(xy[,1]~s(t,bs="cr",k=knots))
  fity <- gam(xy[,2]~s(t,bs="cr",k=knots))

  t2 <- orthogonal_path(fitx,fity,t)

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
                                high="firebrick1")
  p <- p + theme_bw()

  xyt <- data.frame(x=xy[,1],y=xy[,2],t=t)

  p2 <- data.frame(x=xy[,1],y=xy[,2],color=t2) |>
    ggplot(aes(x=x,y=y,color=color)) + geom_point() +
    scale_color_gradient2(low="blue", mid="grey", high="red")

  out <- list()
  out$xyt <- xyt
  out$curve.plot <- p
  out$t2 <- t2
  out$p2 <- p2
  return(out)
}

orthogonal_path <- function(fitx,fity, t) {
  f2x <- gratia::derivatives(fitx,order=2,data=data.frame(t=t))
  f2y <- gratia::derivatives(fity,order=2,data=data.frame(t=t))

  t2 <- t
  first.t <- order(t)[1]
  e <- c(fitx$residuals[first.t], fity$residuals[first.t])
  f2 <- c(f2x$.derivative[first.t], f2y$.derivative[first.t])
  sign <- ifelse(sum(e*f2) > 0, 1, -1)
  t2[first.t] <- sign*sqrt(sum(e^2))
  f2.prev <- f2
  for(i in order(t)[-1]) {
    e <- c(fitx$residuals[i], fity$residuals[i])
    f2 <- c(f2x$.derivative[i], f2y$.derivative[i])

    #if(cosine(f2,f2.prev) > 0.5) {
    # my.sign <- sign(sum(f2*e))
    #} else if(cosine(f2,f2.prev) < -0.5) {
    #  my.sign <- -sign(sum(f2*e))
    #} else{
    #  my.sign <- 0
    #}
    f2 <- sign(sum(f2*f2.prev))*f2
    if(cosine(f2,f2.prev) < 0.3) {
      print("HI")
    }
    my.sign <- sign(sum(f2*e))
    t2[i] <- my.sign*sqrt(sum(e^2))
    f2.prev <- f2
  }

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

