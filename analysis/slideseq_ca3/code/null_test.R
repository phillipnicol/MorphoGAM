
library(STexampleData)
library(dimRed)
library(igraph)
library(tidyverse)
library(splines)
library(mgcv)

source("../../../R/curvesearcher.R")

spe <- SlideSeqV2_mouseHPC()


ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]

#out <- CurveSearcher(xy,knn=5) Use this to fit the curve

Y.sub <- counts(spe)[,ixs]

load("../data/curve_object.rda")


set.seed(1)
#Scramble spatial locations
Y.sub <- Y.sub[,sample(1:ncol(Y.sub), size=ncol(Y.sub), replace=FALSE)]

#First test using GAM
fit <- detectSVG(Y.sub,out)

ngene <- length(which(rowSums(Y.sub) > 20))
qs <- seq(0.001, 0.1, length.out=100)
error_rate <- qs
for(q in 1:length(qs)) {
  error_rate[q] <- sum(fit$p.val < qs[q])/ngene
}

#Now test using ns
fit <- detectSVGns(Y.sub,out)
qs <- seq(0.001, 0.1, length.out=100)
error_rate2 <- qs
for(q in 1:length(qs)) {
  error_rate2[q] <- sum(fit$p.val < qs[q])/ngene
}

df <- data.frame(gam=error_rate, ns=error_rate2)
saveRDS("../data/null_error_rate.RDS")

