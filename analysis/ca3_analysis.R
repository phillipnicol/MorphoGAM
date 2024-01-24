
library(STexampleData)

spe <- SlideSeqV2_mouseHPC()


ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]

out <- CurveSearcher(xy,knn=5)

Y.sub <- counts(spe)[,ixs]

fit1 <- detectSVG(Y.sub, out)


#Scramble spatial locations
Y.sub <- Y.sub[,sample(1:ncol(Y.sub), size=ncol(Y.sub), replace=FALSE)]
fit <- detectSVG(Y.sub,out)

qs <- seq(0.001, 0.1, length.out=100)
error_rate <- qs
for(q in 1:length(qs)) {
  error_rate[q] <- sum(fit$p.val < qs[q])/9225
}



## Simulated example

df <- as.data.frame(xy[-out$outlier,])
sig <- -c(10, 25, 50, 100, 250, 500, 1000)
res <- array(0, dim=c(length(sig), 100, 2))
peak <- res
for(i in 1:length(sig)) {
  alpha <- 0
  true.f <- sig[i]*(out$t - 0.5)^2+1+alpha

  for(j in 1:100) {
    y <- rnegbin(n = length(true.f), mu=n*exp(true.f), theta=5)
    fit <- gam(y~s(out$t)+offset(log(n)),family=nb()) ## TO DO ADD CYCLE
    fx <- fit$linear.predictors - log(n)
    res[i,j,1] <- sqrt(mean(fx - true.f)^2)
    est <- which.max(fx)
    peak[i,j,1] <- abs(out$t[est] - 0.5)

    fit2 <- gam(y~s(df$xcoord, df$ycoord) + offset(log(n)),family=nb())
    fx <- fit2$linear.predictors - log(n)
    res[i,j,2] <- sqrt(mean(fx - true.f)^2)
    est <- which.max(fx)
    peak[i,j,2] <- abs(out$t[est] - 0.5)
  }
}

pow1 <- 0
pow2 <- 0
pow3 <- 0
df <- xy[-out$outlier,]
for(i in 1:20) {
  y <- rnegbin(n = nrow(df), mu=1+1*exp(-1000*(out$t - 0.5)^2), theta=5)
  y <- matrix(y,nrow=1)
  y <- rbind(y,y)
  spark <- sparkx(y, locus_in=df)
  if(spark$res_mtest[1,1] < 0.05/20000) {
    pow1 <- pow1 + 1
  }

  fit <- gam(y[1,]~s(out$t),family=nb()) ## TO DO ADD CYCLE
  if(summary(fit)$s.pv < 0.05/20000) {
    pow2 <- pow2 + 1
  }

  fit3 <- nnSVG(input=log(y+1), spatial_coords = df, verbose=TRUE,
                order="Sum_coords")
  if(fit3[1,13] < 0.05/20000) {
    pow3 <- pow3 + 1
  }
}

## Plotting
df.res <- reshape2::melt(res)
df.peak <- reshape2::melt(peak)

method <- c("1D", "2D")
df.res <- df.res |> group_by(Var1,Var3) |> summarise(val=mean(value))
p <- ggplot(data=df.res, aes(x=-sig[Var1], y=val,color=method[Var3]))
p <- p + geom_point() + geom_line()
p1 <- p + xlab("Sigma") + ylab("MSE of log(mu)")


method <- c("1D", "2D")
df.peak <- df.peak |> group_by(Var1,Var3) |> summarise(val=mean(value))
p <- ggplot(data=df.peak, aes(x=-sig[Var1], y=val,color=method[Var3]))
p <- p + geom_point() + geom_line()
p2 <- p + xlab("Sigma") + ylab("MSE of estimated peak")

library(ggpubr)
p <- ggarrange(p1, p2, nrow=1)


#1000 spots 5 genes
xy <- matrix(runif(2*1000), nrow=1000)
z <- matrix(rnorm(5*1000), nrow=5)

fit <- nnSVG(input=z, spatial_coords = xy, verbose=TRUE)


fit <- nnSVG(input=z, spatial_coords = df, verbose=TRUE, order="Sum_coords")
