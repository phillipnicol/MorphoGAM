load("../data/curve_object.rda")

library(STexampleData)
library(MASS)
library(nnSVG)
library(SPARK)
library(mgcv)

set.seed(1)

spe <- SlideSeqV2_mouseHPC()


ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]

df <- xy[-out$outlier,]

kappa <- c(0.5, 1, 1.5, 2, 2.5, 3)
sigma <- seq(100, 1000, by=100)
niter <- 100

res <- array(0, dim=c(length(kappa), length(sigma), 3))

mu <- 1
theta <- 5
for(i in 1:length(kappa)) {
  for(j in 1:length(sigma)) {
    pow1 <- 0; pow2 <- 0; pow3 <- 0
    for(k in 1:niter) {
      y <- rnegbin(n = nrow(df), mu=mu+kappa[i]*exp(-sigma[j]*(out$t - 0.5)^2), theta=theta)
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
    res[i,j,1] <- pow1/niter
    res[i,j,2] <- pow2/niter
    res[i,j,3] <- pow3/niter
  }
}

saveRDS(res, "../data/results.RDS")

res <- readRDS("../data/results.RDS")

library(reshape2)
library(tidyverse)

df <- melt(res); colnames(df) <- c("kappa", "sigma", "method", "power")

df$kappa <- paste("k = ", kappa[df$kappa])
df$sigma <- sigma[df$sigma]

df$method <- c("SPARK-X", "Projection", "nnSVG")[df$method]

p <- ggplot(data=df,aes(x=sigma, y=power, color=method)) +
  geom_point() +
  geom_line() +
  facet_wrap(~kappa)+
  labs(color="Method")+
  xlab(expression(sigma)) + ylab("Power")+
  theme_bw()
p

ggsave("../plots/sim_power.png")



