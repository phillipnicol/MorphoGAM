

library(STexampleData)
library(igraph)
library(tidyverse)

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

fit <- CurveFinder(xy)

Y <- as.matrix(Y)
mgam <- MorphoGAM(Y, curve.fit=fit,
                  design=y~s(t,bs="cr"))

save(mgam, file="../data/mgam_ca3.RData")
