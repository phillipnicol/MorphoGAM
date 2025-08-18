setwd(here::here("analysis/moffitt_mucosa/code"))

Y <- read.csv("../../data/GL2_distal_colon_cell_by_gene_raw.csv")

rownames(Y) <- Y[,1]

Y <- Y[,-1]

library(Matrix)

Y <- as.matrix(Y)

Y <- t(Y) #Transpose to get genes x cells

meta <- read.csv("../../data/GL2_distal_colon_cell_type_and_locations_2023.08.11.csv")

library(tidyverse)
meta <- as.data.frame(meta)
meta.sub <- meta |> filter(slice_full_name == "20220518_WT_dcol_slice_3") |>
  filter(spatial_neighborhood_v1 == "Mucosa") |>
  filter(leiden_combined_v2 == "Enterocyte")

Y.sub <- Y[,meta.sub$X]

Y.sub <- Y.sub[rowSums(Y.sub) >= 10,]


### SPatial clustering

library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)
library(tidyverse)

counts <- Y.sub |> as.matrix()
colData <- meta.sub
colData <- colData |> mutate(row=x, col=y, array_row=x,array_col=y)
rownames(colData) <- colnames(counts)

sce <- SingleCellExperiment(assays = list(counts = counts),
                             colData = colData)

ngenes <- nrow(Y.sub)

set.seed(102)
sce <- spatialPreprocess(sce, platform="ST",
                              n.PCs=7, n.HVGs=ngenes, log.normalize=TRUE)


sce <- qTune(sce, qs=seq(2, 10), platform="ST")

qPlot(sce)

sce <- spatialCluster(sce, q=3)

clusterPlot(sce)
