
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
fit1 <- detectSVG(Y.sub,out)

#Now test using ns
fit2 <- detectSVGns(Y.sub,out)

df <- data.frame(gam=fit1$p.val, ns=fit2$p.val)
saveRDS("../data/null_error_rate.RDS")

