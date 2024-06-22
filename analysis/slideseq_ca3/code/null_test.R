
library(STexampleData)
library(tidyverse)
library(splines)
library(mgcv)


spe <- STexampleData::SlideSeqV2_mouseHPC()

ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]
Y <- counts(spe)[,ixs]

xy.dist <- as.matrix(dist(xy))
knn <- 20
prune.outlier <- 3
nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
outlier <- which(nnk > (prune.outlier)*median(nnk))

xy <- xy[-outlier,]; Y <- Y[,-outlier]

load("../data/ca3_curve.RData")
t <- fit$xyt$t

set.seed(1)
#Scramble spatial locations
Y <- Y[,sample(1:ncol(Y), size=ncol(Y), replace=FALSE)]
Y <- Y[rowSums(Y) >= 10,]
l.o <- log(colSums(Y))

p.val <- rep(0, nrow(Y))
#First test using GAM
for(j in 1:nrow(Y)) {
  print(j)
  fit1 <- mgcv::gam(Y[j,] ~ s(t,bs="cr") + offset(l.o),family=nb(), H=diag(10))
  p.val[j] <- summary(fit1)$s.pv
}

saveRDS(p.val, file= "../data/null_error_rate.RDS")


## Type I error
cutoff <- seq(10^{-4}, 0.1, length.out=1000)
t1e <- sapply(cutoff, function(i) sum(p.val < i))/nrow(Y)

p <- data.frame(x=cutoff, y=t1e) |>
  ggplot(aes(x=x,y=y)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
  xlab("Significance level") +
  ylab("Type I error rate") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

ggsave(filename="../plots/type_1_error.png", width=6.1, height=4.85, units="in")

