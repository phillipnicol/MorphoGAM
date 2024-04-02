
##Integral of bias^2

#Mean of square bias

## Integral of variance
#Approximately sigma^2/nbhd size

CurveSearcher.cv <- function(xy) {
  #Try
  knn <- c(3,5,10,15,20)
  tau <- c(10,50,100,250)
  cutoff <- c(2,5,10)
  model <- c(1,2)

  params <- expand.grid(knn,tau,cutoff,model)
  params <- as.matrix(params)
  colnames(params) <- c("knn", "tau",
                        "cutoff", "model")
  mse <- apply(params, 1, function(param) {
    if(param["model"] == 1) {
      fit <- CurveSearcher(xy,
                           knn=param["knn"],
                           tau=param["tau"],
                           cutoff=param["cutoff"])
      fit$mse
    } else {
      fit <- CurveSearcherLoop(xy,
                           knn=param["knn"],
                           tau=param["tau"],
                           cutoff=param["cutoff"])
      fit$mse
    }
  })
}

bias <- function(xy, ft) {
  error <- rep(0, nrow(xy))
  for(j in 1:nrow(xy)) {
    my.dist <- apply(out$ft, 1, function(x) sum((xy[j,] - x)^2))
    predicted.xy <- out$ft[which.min(my.dist),]
    error[j] <- sum((xy[j,] - predicted.xy)^2)
  }
  mean(error)

  sigma.hat <- sd(error)
  avg.var <- 2*sigma.hat^2*mean(out$my.var)
  avg.var

  mean(error)+avg.var
}




















CurveSearcher <- function(xy, knn, tau=100,
                          cutoff=2) {

  xy.dist <- as.matrix(dist(xy))

  nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
  outlier <- which(nnk > cutoff*median(nnk))

  if(length(outlier) > 0) {
    xy.new <- xy[-outlier,]
  } else {
    xy.new <- xy
  }

  knng <- dimRed:::makeKNNgraph(x = xy.new,
                                k = knn,
                                eps = 0)


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
  t <- rep(0, nrow(xy.new)); names(t) <- rownames(xy.new)
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
  my.var <- my.t
  for(i in 1:nrow(ft)) {
    #print(i)

    my.dist <- abs(t - my.t[i])
    kw <- exp(-tau*my.dist)
    ft[i,1] <- weighted.mean(x=xy.new[,1], w=kw)
    ft[i,2] <- weighted.mean(x=xy.new[,2], w=kw)

    my.var[i] <- sum(kw^2)/(sum(kw)^2)
  }

  ## Project the outliers back
  t.new <- rep(0, nrow(xy))
  for(i in 1:nrow(xy)) {
    if(rownames(xy)[i] %in% names(t)) {
      ix <- which(rownames(xy)[i] == names(t))
      t.new[i] <- t[ix]
    } else {
      my.dist <- apply(ft, 1, function(x) sum((xy[i,] - x)^2))
      t.new[i] <- my.t[which.min(my.dist)]
    }
  }

  error <- rep(0, nrow(xy))
  error.x <- error
  error.y <- error
  for(j in 1:nrow(xy)) {
    my.dist <- apply(ft, 1, function(x) sum((xy[j,] - x)^2))
    predicted.xy <- ft[which.min(my.dist),]
    error[j] <- sum((xy[j,] - predicted.xy)^2)
    error.x[j] <- (xy[j,1] - predicted.xy[1])^2
    error.y[j] <- (xy[j,2] - predicted.xy[2])^2
  }

  df <- sum(my.var)
  sigma.x <- sum(error.x)/(nrow(xy)-df)
  sigma.y <- sum(error.y)/(nrow(xy)-df)
  avg.var <- (sigma.x+sigma.y)*mean(my.var)


  null.bias <- (xy[,1] - mean.x)^2 + (xy[,2] - mean.y)^2
  null.bias <- mean(null.bias)
  sigma.hat.null <- 0.5*(var(xy[,1]) + var(xy[,2]))
  null.var <- 2*sigma.hat.null/nrow(xy)
  null.mse <- null.bias + null.var


  df.new <- data.frame(x=xy[,1],
                       y=xy[,2])

  p <- ggplot(data=df.new,aes(x=x,y=y)) + geom_point(col="grey",
                                                     alpha=0.5)

  df.line <- data.frame(x=ft[,1], y=ft[,2], color=seq(0,1,by=0.001))
  p <- p + geom_path(data=df.line,aes(x=x,y=y,color=color),
                     linewidth=1)
  p <- p + scale_color_gradient(low="navyblue",
                                high="firebrick1")
  p <- p + theme_bw()

  xy.dist <- as.matrix(dist(xy))
  t.dist <- as.matrix(dist(t.new))

  out <- list()
  out$plot <- p
  out$t <- t
  out$df <- df
  out$t.new <- t.new
  out$outlier <- outlier
  out$my.var <- my.var
  out$ft <- ft
  out$R.hat <- 1 - (mean(error)+avg.var)/null.mse
  out$null.mse <- null.mse
  out$mse <- mean(error)+avg.var
  return(out)
}







