load("../data/ca3_curve.RData")
out <- fit

library(STexampleData)
library(MASS)
library(nnSVG)
library(SPARK)
library(mgcv)
library(splines)

set.seed(1)


spe <- STexampleData::SlideSeqV2_mouseHPC()

ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]
Y <- counts(spe)[,ixs]

xy.dist <- as.matrix(dist(xy))
knn <- 20
prune.outlier <- 3
nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
outlier <- which(nnk > (prune.outlier)*median(nnk))

df <- xy[-outlier,]




kappa <- seq(0.5, 1.5, length.out=6)
sigma <- seq(100, 1000, by=100)
niter <- 100

res <- array(0, dim=c(length(kappa), length(sigma),4))

mu <- 1
theta <- 5
for(i in 1:length(kappa)) {
  for(j in 1:length(sigma)) {
    print(i); print(j); 
    pow1 <- 0; pow2 <- 0; pow3 <- 0; pow4 <- 0
    for(k in 1:niter) {
      y <- rnegbin(n = nrow(df), mu=(exp(-sigma[j]*(out$xyt$t - 0.5)^2))*kappa[i] + 1,
                   theta=theta)
      #y <- rnegbin(n = nrow(df), mu=2*sin(40*out$t)+2, theta=theta)
      y <- matrix(y,nrow=1)
      y <- rbind(y,y)
      spark <- sparkx(y, locus_in=df)
      if(spark$res_mtest[1,1] < 0.05/20000) {
        pow1 <- pow1 + 1
      }

      
      fit <- gam(y[1,]~s(out$xyt$t,k=10,bs="cr"),family=nb(),
                 H=diag(10))
      if(summary(fit)$s.pv < 0.05/20000) {
        pow2 <- pow2 + 1
      }

      fit3 <- nnSVG(input=log(y+1), spatial_coords = df, verbose=TRUE,
                    order="Sum_coords")
      if(fit3[1,13] < 0.05/20000) {
        pow3 <- pow3 + 1
      }

      #fit4 <- glm(y[1,] ~ ns(out$t,df=10), family=quasipoisson())
      #fx <- fit4$linear.predictors - fit4$coefficients[1]
      #fx <- fx - mean(fx)
      #Sigma <- vcov(fit4)[-1,-1]
      #Z <- t(rmvnorm(n=10^4,sigma=Sigma))
      #X <- model.matrix(fit4)
      #fx.null <- X[,-1] %*% Z
      #fx.null <- sweep(fx.null, MARGIN=2,
      #                 STATS=colMeans(fx.null), FUN = "-")
      #peaks <- apply(fx.null, 2, max)
      #p.val <- (sum(peaks > max(fx)) + 1)/(10^4 + 1)

      #if(p.val < 0.05/20000) {
      #  pow4 <- pow4+1
      #}
    }
    res[i,j,1] <- pow1/niter
    res[i,j,2] <- pow2/niter
    res[i,j,3] <- pow3/niter
    res[i,j,4] <- pow4/niter
  }
}

saveRDS(res, "../data/results.RDS")

res <- readRDS("../data/results.RDS")
res <- res[,,-4]
library(reshape2)
library(tidyverse)
kappa <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
sigma <- seq(100, 1000, by=100)
df <- melt(res); colnames(df) <- c("kappa", "sigma", "method", "power")
df$kappa <- paste("k = ", kappa[df$kappa])
df$sigma <- sigma[df$sigma]

df$method <- c("SPARK-X", "Projection", "nnSVG")[df$method]

p.power <- ggplot(data=df,aes(x=sigma, y=power, color=method)) +
  geom_point() +
  geom_line() +
  facet_wrap(~kappa)+
  labs(color="Method")+
  xlab(expression(sigma)) + ylab("Power")+
  theme_bw() + theme(legend.position = "top")


## Prune outlier


library(STexampleData)
spe <- SlideSeqV2_mouseHPC()
ixs <- which(spe$celltype == "CA3") #subset to CA3
xy <- spatialCoords(spe)[ixs,]

xy.dist <- as.matrix(dist(xy))
nnk <- apply(xy.dist, 1, function(x) sort(x)[10])
outlier <- which(nnk > 2*median(nnk))
xy <- xy[-outlier,]

fit <- CurveFinder(xy)

p.curve <- fit$curve.plot + guides(color="none")+
  xlab("") + ylab("") + ggtitle("Fitted curve")

p.coord <- fit$coordinate.plot + guides(color="none") +
  xlab("") + ylab("")

p.resid <- fit$residuals.plot + guides(color="none") +
  xlab("") + ylab("")

p <- ggarrange(ggarrange(p.curve,p.coord,p.resid,nrow=1),
          p.power, nrow=2,ncol=1,
          heights=c(1,1.5))

ggsave(p,filename="../plots/sim_power.png",
       width=7,height=10)



