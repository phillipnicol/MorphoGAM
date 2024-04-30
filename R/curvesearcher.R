

PathFinder <- function(xy) {
  #Try
  knn <- c(3,5,10,15,20)
  cutoff <- c(2,5,10)
  model <- c(1,2)

  params <- expand.grid(knn,cutoff,model)
  params <- as.matrix(params)
  colnames(params) <- c("knn","cutoff", "model")
  mse <- apply(params, 1, function(param) {
    if(param["model"] == 1) {
      print(param)
      fit <- CurveSearcher(xy,
                           knn=param["knn"],
                           cutoff=param["cutoff"])
      fit$mse
    } else {
      print(param)
      fit <- CurveSearcherLoop(xy,
                               knn=param["knn"],
                               cutoff=param["cutoff"])
      fit$mse
    }
  })
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

  knng <- dimRed:::makeKNNgraph(x = xy,
                                k = knn,
                                eps = 0)


  comp <- components(knng)


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


  best.path <- which.min(paths$cost)
  t.shift <- 0
  t <- rep(0, nrow(xy))
  for(c in 1:comp$no) {
    current.comp <- paths[best.path,c]

    if(current.comp > 0) {
      t[comp$membership == current.comp] <- t.list[[current.comp]] + t.shift
    } else{
      t.max <- max(t.list[[-current.comp]])
      t[comp$membership == -current.comp] <- t.max - t.list[[-current.comp]] + t.shift
    }
    t.shift <- max(t)
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


  knots <- min(100, 0.1*nrow(xy.new))
  fitx <- gam(xy.new[,1]~s(t,bs="cr",k=knots))
  fity <- gam(xy.new[,2]~s(t,bs="cr",k=knots))

  ## Project the outliers back
  t.new <- rep(0, nrow(xy))
  for(i in 1:nrow(xy)) {
    if(rownames(xy)[i] %in% names(t)) {
      ix <- which(rownames(xy)[i] == names(t))
      t.new[i] <- t[ix]
    } else {
      knnt <- order(xy.dist[i,-outlier], decreasing=FALSE)[1:knn]
      t.new[i] <- mean(t[knnt])
    }
  }

  my.t <- seq(0,1,by=10^{-4})
  predx <- predict(fitx,newdata=list(t=my.t))
  predy <- predict(fity,newdata=list(t=my.t))

  fp <- diff(predx)
  fp2 <- diff(predy)
  al <- sqrt(fp^2 + fp2^2)

  ft <- data.frame(x=predx,y=predy)


  SStot <- sum((xy[,1]-mean(xy[,1]))^2+(xy[,2]-mean(xy[,2]))^2)
  SSresid <- sum(fitx$residuals^2 + fity$residuals^2)

  ft.fitted <- cbind(predict(fitx,newdata=list(t=t.new)),
                 predict(fity,newdata=list(t=t.new)))
  resid <- (xy[,1] - ft.fitted[,1])^2 + (xy[,2] - ft.fitted[,2])^2
  bias <- median(resid)

  # Now compute variance

#  ft.fitted <- cbind(predict(fitx, newdata=list(t=t.new)),
#                     predict(fity, newdata=list(t=t.new)))
#  my.var <- rep(0, nrow(xy))
#  for(i in 1:nrow(xy)) {
#    ft.i <- ft.fitted[i,]
#    knn.i <- order(xy.dist[i,], decreasing=FALSE)[1:knn]
#    varx <- var(ft.fitted[knn.i,1])
#    vary <- var(ft.fitted[knn.i,2])
#    my.var[i] <- varx + vary
#  }

  my.var <- rep(0, length(my.t))
  for(i in 1:length(my.t)) {
    range <- which(abs(my.t-my.t[i]) < 0.01)
    my.var[i] <- var(ft[range,1]) + var(ft[range,2])
  }



  #ft <- matrix(0, nrow=length(my.t), ncol=2)
  #my.var <- my.t
  #for(i in 1:nrow(ft)) {
    #print(i)

  #  my.dist <- abs(t - my.t[i])
  #  kw <- exp(-tau*my.dist)
  #  ft[i,1] <- weighted.mean(x=xy.new[,1], w=kw)
  #  ft[i,2] <- weighted.mean(x=xy.new[,2], w=kw)

  #  my.var[i] <- sum(kw^2)/(sum(kw)^2)
  #}


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

  out <- list()
  out$plot <- p
  out$t <- t.new
  out$outlier <- outlier
  out$ft <- ft
  out$al <- sum(al)
  out$Rhat <- 1-SSresid/SStot
  out$fitx <- fitx; out$fity <- fity
  out$mse <- bias + mean(my.var)
  out$bias <- bias
  out$my.var <- my.var
  return(out)
}



CurveSearcherLoop <- function(xy, knn, tau=100,
                              cutoff=2) {

  xy.dist <- as.matrix(dist(xy))

  nn1 <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
  outlier <- which(nn1 > cutoff*median(nn1))

  if(length(outlier) > 0) {
    xy.new <- xy[-outlier,]
  } else {
    xy.new <- xy
  }

  knng <- dimRed:::makeKNNgraph(x = xy.new,
                                k = knn,
                                eps = 0)

  comp <- components(knng)
  if(comp$no > 1) {
    out <- list()
    out$mse <- Inf
    out$Rhat <- 0
    out$al <- 1
    return(out)
  }

  geodist <- igraph::distances(knng, algorithm = "dijkstra")

  k <- geodist ^ 2
  k <- .Call(stats:::C_DoubleCentre, k)
  k <- - k / 2

  e <- RSpectra::eigs_sym(k, 4, which = "LA",
                          opts = list(retvec = TRUE))

  t <- (pi+atan2(e$vectors[,1], e$vectors[,2]))/(2*pi)
  r <- mean(sqrt(e$values[1]*e$vectors[,1]^2 + e$values[2]*e$vectors[,2]^2))
  names(t) <- rownames(xy.new)

  my.t <- t
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

  knots <- min(100, 0.1*nrow(xy.new))
  fitx <- gam(xy.new[,1]~s(t,bs="cc",k=knots))
  fity <- gam(xy.new[,2]~s(t,bs="cc",k=knots))

  ## Project the outliers back
  t.new <- rep(0, nrow(xy))
  for(i in 1:nrow(xy)) {
    if(rownames(xy)[i] %in% names(t)) {
      ix <- which(rownames(xy)[i] == names(t))
      t.new[i] <- t[ix]
    } else {
      knnt <- order(xy.dist[i,-outlier])[1:knn]
      t.new[i] <- mean(t[knnt])
    }
  }

  my.t <- seq(0,1,by=10^{-4})
  predx <- predict(fitx,newdata=list(t=my.t))
  predy <- predict(fity,newdata=list(t=my.t))


  fp <- diff(predx)
  fp2 <- diff(predy)
  al <- sqrt(fp^2 + fp2^2)

  SStot <- sum((xy[,1]-mean(xy[,1]))^2+(xy[,2]-mean(xy[,2]))^2)
  SSresid <- sum(fitx$residuals^2 + fity$residuals^2)

  ft.fitted <- cbind(predict(fitx,newdata=list(t=t.new)),
                 predict(fity,newdata=list(t=t.new)))
  resid <- (xy[,1] - ft.fitted[,1])^2 + (xy[,2] - ft.fitted[,2])^2

  bias <- mean(resid)

  ft <- data.frame(x=predx,y=predy)


  my.var <- rep(0, length(my.t))
  for(i in 1:length(my.t)) {
    range <- which(abs(my.t-my.t[i]) < 0.01)
    my.var[i] <- var(ft[range,1]) + var(ft[range,2])
  }

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

  xy.dist <- as.matrix(dist(xy))
  t.dist <- as.matrix(dist(t))
  out <- list()
  out$plot <- p
  out$t <- t.new
  out$outlier <- outlier
  out$ft <- ft
  out$r <- r
  out$al <- sum(al)
  out$Rhat <- 1-SSresid/SStot
  out$fitx <- fitx; out$fity <- fity
  out$bias <- bias
  out$var <- my.var
  out$mse <- bias+mean(my.var)
  return(out)
}

CurveSearcher.cv <- function(xy, nfolds=5) {
  my.folds <- createFolds(1:nrow(xy), k = nfolds)

  knn <- c(3,5,10,15,20)
  cutoff <- c(2,5,10)
  model <- c(1,2)

  params <- expand.grid(knn,cutoff,model)
  params <- as.matrix(params)
  colnames(params) <- c("knn","cutoff", "model")

  mse <- matrix(0, nrow=nfolds, ncol=nrow(params))

  for(k in 1:nfolds) {
    xy.sub <- xy[-my.folds[[k]],]
    xy.ho <- xy[my.folds[[k]],]

    mse.k <- apply(params, 1, function(param) {
      if(param["model"] == 1) {
        print(param)
        fit <- CurveSearcher(xy.sub,
                             knn=param["knn"],
                             cutoff=param["cutoff"])
      } else {
        print(param)
        fit <- CurveSearcherLoop(xy.sub,
                                 knn=param["knn"],
                                 cutoff=param["cutoff"])
      }

      print(fit$Rhat)
      if(fit$Rhat < 10^{-10}) {
        return(Inf)
      }
      error <- rep(0, nrow(xy.ho))
      for(j in 1:nrow(xy.ho)) {
        my.dist <- apply(xy.sub, 1, function(x) sum((xy.ho[j,] - x)^2))
        predicted.t <- mean(fit$t[order(my.dist)[1:param["knn"]]])
        pred.x <- predict(fit$fitx, newdata=list(t=predicted.t))
        pred.y <- predict(fit$fity, newdata=list(t=predicted.t))
        predicted.xy <-  c(pred.x, pred.y)
        error[j] <- sum((xy.ho[j,] - predicted.xy)^2)
      }
      return(sqrt(mean(error)))
    })

    mse[k,] <- mse.k
  }

  mse.avg <- colMeans(mse)
  best.param <- params[which.min(mse.avg),]

  if(best.param["model"] == 1) {
    fit <- CurveSearcher(xy,
                         knn=best.param["knn"],
                         cutoff=best.param["cutoff"])
  } else{
    fit <- CurveSearcherLoop(xy,
                             knn=best.param["knn"],
                             cutoff=best.param["cutoff"])
  }
  return(fit)
}


CurveSearcher.cv2 <- function(xy) {
  xy.sub <- xy
  hold.out <- sample(1:nrow(xy.sub), size=300,replace=FALSE)
  xy.ho <- xy.sub[hold.out,]
  xy.sub <- xy.sub[-hold.out,]

  knn <-10
  out <- CurveSearcher(xy.sub, knn=knn,cutoff=2)


  knn <- 20
  out <- CurveSearcherLoop(xy.sub,knn=knn,cutoff=10)
  error <- rep(0, nrow(xy.ho))
  for(j in 1:nrow(xy.ho)) {
    my.dist <- apply(xy.sub, 1, function(x) sum((xy.ho[j,] - x)^2))
    predicted.t <- mean(out$t[order(my.dist)[1:knn]])
    pred.x <- predict(out$fitx, newdata=list(t.new=predicted.t))
    pred.y <- predict(out$fity, newdata=list(t.new=predicted.t))
    predicted.xy <-  c(pred.x, pred.y)
    error[j] <- sum((xy.ho[j,] - predicted.xy)^2)
  }
  sqrt(mean(error))
  quantile(sqrt(error),0.75)

  tau_lb <- 0; tau_ub <- 10^3

  tau.try <- seq(10, 1000, by=50)
  rmse <- rep(0, length(tau.try))
  for(i in 1:length(tau.try)) {
    out <- CurveSearcher(xy.sub,knn=2, tau=tau.try[i])

    error <- rep(0, nrow(xy.ho))
    for(j in 1:nrow(xy.ho)) {
      my.dist <- apply(out$ft, 1, function(x) sum((xy.ho[j,] - x)^2))
      predicted.xy <- out$ft[which.min(my.dist),]
      error[j] <- sum((xy.ho[j,] - predicted.xy)^2)
    }

    rmse[i] <- sqrt(median(error))
    cat(tau.try[i], " ", rmse[i], "\n")
  }

  error <- rep(0, nrow(xy.ho))
  for(j in 1:nrow(xy.ho)) {
    my.dist <- apply(out$ft, 1, function(x) sum((xy.ho[j,] - x)^2))
    predicted.xy <- out$ft[which.min(my.dist),]
    error[j] <- sum((xy.ho[j,] - predicted.xy)^2)
  }

  lb_rmse <- sqrt(mean(error))

  out <- CurveSearcher(xy.sub,knn=knn, tau=10^3)

  error <- rep(0, nrow(xy.ho))
  for(j in 1:nrow(xy.ho)) {
    my.dist <- apply(out$ft, 1, function(x) sum(xy.ho[1,] - x)^2)
    predicted.xy <- out$ft[which.min(my.dist),]
    error[j] <- sum((xy.ho[j,] - predicted.xy)^2)
  }

  ub_rmse <- sqrt(mean(error))

  for(knn in 5:15) {
    while(TRUE) {
      out <- CurveSearcher(xy.sub,knn=knn, tau=tau)
      error <- rep(0, nrow(xy.ho))
      for(j in 1:nrow(xy.ho)) {
        my.dist <- apply(out$ft, 1, function(x) sum(xy.ho[1,] - x)^2)
        predicted.xy <- out$ft[which.min(my.dist),]
        error[j] <- sum((xy.ho[j,] - predicted.xy)^2)
      }

      rmse <- sqrt(mean(error))
      cat(tau, " ", rmse, "\n")
      if(rmse > tau.rmse.prev & (tau-tau.start)/tau < 0.1) {
        break
      } else if(rmse > tau.rmse.prev) {
        old.tau <- tau
        tau <- (tau+tau.start)/2
        tau.start <- old.tau
      } else {
        tau <- 2*tau
        tau.start <- 1/2*tau
      }

      tau.rmse.prev <- rmse
    }
    out <- CurveSearcher(xy.sub,knn=knn, tau=tau.start)

    error <- rep(0, nrow(xy.ho))
    for(j in 1:nrow(xy.ho)) {
      my.dist <- apply(out$ft, 1, function(x) sum(xy.ho[1,] - x)^2)
      predicted.xy <- out$ft[which.min(my.dist),]
      error[j] <- sum((xy.ho[j,] - predicted.xy)^2)
    }


  }
}

kernelSmoother <- function(xy, t, tau) {
  my.t <- seq(0, 1, by=0.001)
  ft <- matrix(0, nrow=length(my.t), ncol=2)
  for(i in 1:nrow(ft)) {
    #print(i)

    my.dist <- abs(t - my.t[i])
    kw <- exp(-tau*my.dist)
    ft[i,1] <- weighted.mean(x=xy.sub[,1], w=kw)
    ft[i,2] <- weighted.mean(x=xy.sub[,2], w=kw)
  }
  return(ft)
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

    fit <- glm.nb(gene~ns(t,df=10)+offset(l.o))
    fit <- update(fit, method = "brglmFit", type = "correction")

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

createPenalty <- function(X.model, cutoff=log(2)) {
  nsim <- 10^4
  p <- ncol(X.model)
  Z <- matrix(rnorm(nsim*p),nrow=p, ncol=nsim)
  fx.null <- X.model %*% Z
  peaks <- apply(fx.null, 2, max)
  q97.5 <- quantile(peaks,0.975)
  sigma <- log(2)/q97.5
  return(1/(2*sigma^2))
}

detectSVGLoop <- function(Y, cso) {
  t <- cso$t; outlier <- cso$outlier
  l.o <- log(colSums(Y[,-outlier]))

  res <- data.frame(peak = rep(0,nrow(Y)),
                    range = rep(0,nrow(Y)),
                    p.val = rep(1, nrow(Y)))

  gene <- rep(0, length(t))
  G <- gam(gene~s(t,k=10,bs="cr"), family=nb(), fit=FALSE)
  lambda <- createPenalty(G$X[,-1])
  f <- matrix(0, nrow=nrow(Y), ncol=ncol(Y))
  for(k in 1:nrow(res)) {
    gene <- Y[k,]

    H.pen <- lambda*diag(10); H.pen[1,1] <- 0
    fit <- gam(gene~s(t,k=10, bs="cr")+offset(l.o),family=nb(),
               H=H.pen)

    #fit <- glm(gene~ns(t,df=10)+offset(l.o),family=quasipoisson())

    fx <- fit$linear.predictors - l.o - fit$coefficients[1]
    f[k,] <- fx

    res$peak[k] <- max(fx)
    res$range[k] <- max(exp(fx+fit$coefficients[1])) - min(exp(fx+fit$coefficients[1]))
    res$p.val[k] <- summary(fit)$s.pv

    cat(k, " ", res$peak[k], " ", res$p.val[k], "\n")
  }
  out <- list()
  out$res <- res
  out$f <- f
  return(out)
}



detectSVGns <- function(Y, cso) {
  t <- cso$t; outlier <- cso$outlier
  l.o <- log(colSums(Y[,-outlier]))

  res <- data.frame(peak = rep(0,nrow(Y)),
                    range = rep(0,nrow(Y)),
                    p.val = rep(1, nrow(Y)))

  f <- matrix(0, nrow=nrow(Y), ncol=ncol(Y)-length(outlier))
  for(k in 1:nrow(res)) {
    gene <- Y[k,-outlier]

    if(sum(gene != 0) < 50) {
      next
    }
    fit <- glm(gene ~ ns(cso$t, df=30)+offset(l.o), family=quasipoisson())
    fit.null <- glm(gene ~ 1, family=quasipoisson())
    val <- anova(fit, fit.null, test="LRT")
    p.val <- val$`Pr(>Chi)`[2]

    fx <- fit$linear.predictors - l.o - fit$coefficients[1]
    fx <- fx - mean(fx)
    Sigma <- vcov(fit)[-1,-1]
    Z <- t(rmvnorm(n=10^4,sigma=Sigma))
    fx.null <- X[,-1] %*% Z
    fx.null <- sweep(fx.null, MARGIN=2,
                     STATS=colMeans(fx.null), FUN = "-")

    peaks <- apply(fx.null, 2, max)
    res$p.val[k] <- (sum(peaks > max(fx)) + 1)/(10^4 + 1)

    res$peak[k] <- max(fx)

    cat(k, " ", res$peak[k], " ", res$p.val[k],  "\n")
  }
  out <- list()
  out$res <- res
  out$f <- f
  return(out)
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
          cost <- cost + my.dist[-2*x[j], -2*x[j+1]]
        } else{
          cost <- cost + my.dist[-2*x[j], 2*x[j+1]-1]
        }
      } else{
        if(x[j+1] < 0) {
          cost <- cost + my.dist[2*x[j]-1, -2*x[j+1]]
        } else{
          cost <- cost + my.dist[2*x[j]-1, 2*x[j+1]-1]
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

