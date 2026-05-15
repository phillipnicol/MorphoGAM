setwd(here::here("analysis/moffitt_mucosa/code"))

library(RColorBrewer)
library(ggrepel)
library(ggpubr)

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

old.rownames <- rownames(Y.sub)



for(k in 7:20) {
  fit <- CurveFinder(locus, loop=TRUE,knn=k)
  print(fit$curve.plot)
  #Save model score
  model.score <- fit$model.score
  print(paste("k=",k,"model score=",model.score))
}




for(k in 4:15) {
  fit <- CurveFinder(xy, loop=FALSE,knn=k)
  print(fit$curve.plot)
  #Save model score
  model.score <- fit$model.score
  print(paste("k=",k,"model score=",model.score))
}


  fit <- CurveFinder(xy, loop=FALSE,knn=9)
