


##Make curve
library(tidyverse)

ad <- anndata::read_h5ad("/Users/phillipnicol/Downloads/adata_swiss_2023_10_31_annotated_cleaned.h5ad")

xy.061923 <- ad$obs[ad$obs$Mouse_ID == "061923_D9_m2_Swiss",] |>
  select(x,y)

library(MorphoGAM)

xy <- xy.061923 |> as.matrix()
fit <- CurveFinderInteractive(xy)

save(fit,file="../data/curve_d9_061923.RData")


Y <- ad$raw$X[ad$obs$Mouse_ID == "061923_D9_m2_Swiss",]

library(Matrix)

Y <- as(Y,"sparseMatrix")

Matrix::writeMM(Y,file="../data/d9_061923.mtx")






xy.080823_m5 <- ad$obs[ad$obs$Mouse_ID == "080823_D9_m5_Swiss",] |>
  select(x,y)

library(MorphoGAM)

xy <- xy.080823_m5 |> as.matrix()
fit <- CurveFinderInteractive(xy)

save(fit,file="../data/curve_d9_m5_080823.RData")


Y <- ad$raw$X[ad$obs$Mouse_ID == "080823_D9_m5_Swiss",]

library(Matrix)

Y <- as(Y,"sparseMatrix")

Matrix::writeMM(Y,file="../data/080823_d9_m5.mtx")







xy.080823_m13 <- ad$obs[ad$obs$Mouse_ID == "080823_D9_m13_Swiss",] |>
  select(x,y)

library(MorphoGAM)

xy <- xy.080823_m13 |> as.matrix()
fit <- CurveFinderInteractive(xy)

save(fit,file="../data/curve_d9_m13_080823.RData")


Y <- ad$raw$X[ad$obs$Mouse_ID == "080823_D9_m13_Swiss",]

library(Matrix)

Y <- as(Y,"sparseMatrix")

Matrix::writeMM(Y,file="../data/080823_d9_m13.mtx")


