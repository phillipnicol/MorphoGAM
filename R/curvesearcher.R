

CurveSearcher <- function(xy, knn) {

  xy.dist <- as.matrix(dist(xy))


  nn1 <- apply(xy.dist, 1, function(x) sort(x)[10])
  outlier <- which(nn1 > 3*median(nn1))

  if(length(outlier) > 0) {
    xy.new <- xy[-outlier,]
  } else {
    xy.new <- xy
  }

  #knng <- dimRed:::makeKNNgraph(x = xy.new,
  #                              k = knn,
  #                              eps = 0)

  knng <- makeSNNGraph(xy.new,k=knn)


  comp <- components(knng)

  endpoints_low <- matrix(0, nrow=comp$no, ncol=2)
  rownames(endpoints_low) <- 1:comp$no
  endpoints_high <- matrix(0, nrow=comp$no, ncol=2)

  endpoints <- matrix(0, nrow=2*comp$no, ncol=2)
  t.list <- list()
  for(c in 1:comp$no) {
    xy.sub <- xy.new[comp$membership == c,]

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

    endpoints_low[c,] <- xy.sub[which.min(t.list[[c]]),]
    endpoints_high[c,] <- xy.sub[which.max(t.list[[c]]),]
    endpoints[2*c-1,] <- xy.sub[which.max(t.list[[c]]),]
    endpoints[2*c,] <- xy.sub[which.min(t.list[[c]]),]
  }


  my.dist <- as.matrix(dist(endpoints))
  start_point <- which.max(-endpoints[,1] + endpoints[,2])
  t.shift <- 0
  t <- rep(0, nrow(xy.new))
  for(c in 1:comp$no) {
    current.comp <- ceiling(start_point/2)

    my.dist[,start_point] <- Inf
    if((start_point %% 2) == 0) {
      t[comp$membership == current.comp] <- t.list[[current.comp]] + t.shift
      my.dist[,start_point-1] <- Inf
    } else{
      t.max <- max(t.list[[current.comp]])
      t[comp$membership == current.comp] <- t.max - t.list[[current.comp]] + t.shift
      my.dist[,start_point+1] <- Inf
    }
    t.shift <- max(t)

    if(c == comp$no) {
      break
    }
    start_point <- which.min(my.dist[start_point,])
  }


  my.t <- seq(0, 1, by=0.001)
  ft <- matrix(0, nrow=length(my.t), ncol=2)
  for(i in 1:nrow(ft)) {
    #print(i)

    my.dist <- abs(t - my.t[i])
    kw <- exp(-150*my.dist)
    ft[i,1] <- weightedMedian(x=xy.new[,1], w=kw)
    ft[i,2] <- weightedMedian(x=xy.new[,2], w=kw)
  }

  df.new <- data.frame(x=xy.new[,1],
                       y=xy.new[,2])

  p <- ggplot(data=df.new,aes(x=x,y=y)) + geom_point(col="grey",
                                                     alpha=0.5)

  df.line <- data.frame(x=ft[,1], y=ft[,2], color=seq(0,1,by=0.001))
  p <- p + geom_path(data=df.line,aes(x=x,y=y,color=color),
                     linewidth=1)
  p <- p + scale_color_gradient(low="navyblue",
                                high="firebrick1")
  p <- p + theme_bw()

  out <- list()
  out$plot <- p
  out$t <- t
  out$outlier <- outlier
  return(out)
}

detectSVG <- function(Y, cso) {
  t <- cso$t; outlier <- cso$outlier
  l.o <- log(colSums(Y[,-outlier]))

  res <- data.frame(peak = rep(0,nrow(Y)),
                    range = rep(0,nrow(Y)),
                    p.val = rep(1, nrow(Y)))

  f <- matrix(0, nrow=nrow(Y), ncol=ncol(Y)-length(outlier))
  for(k in 1:nrow(res)) {
    gene <- Y[k,-outlier]

    if(sum(gene != 0) < 20) {
      next
    }

    fit <- gam(gene~s(t)+offset(l.o),family=nb()) ## TO DO ADD CYCLE
    fx <- fit$linear.predictors - l.o - fit$coefficients[1]
    f[k,] <- fx

    res$peak[k] <- max(abs(fx))
    res$range[k] <- max(exp(fx+fit$coefficients[1])) - min(exp(fx+fit$coefficients[1]))
    res$p.val[k] <- summary(fit)$s.pv

    cat(k, " ", res$peak[k], " ", res$p.val[k], " ", res$range[k], "\n")
  }
  out <- list()
  out$res <- res
  out$f <- f
  return(out)
}
