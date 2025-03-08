
setwd(here::here("analysis/merfish_swissroll/code"))

library(anndata)

ad <- anndata::read_h5ad("../../data/adata_swiss_2023_10_31_annotated_clean.h5ad")


meta <- ad$obs

library(MorphoGAM)
library(tidyverse)

meta.sub1 <- meta |> filter(Mouse_ID == "061923_D9_m2_Swiss") |>
  filter(Tier1 %in% c("Fibroblast", "Epithelial"))

xy <- meta.sub1 |> select(x,y) |> as.matrix()
fit <- CurveFinderInteractive(xy)

save(fit, file="../data/061923_D9_m2_Swiss_fibroblast.RData")

# Save matrix

Y <- ad$raw[rownames(meta.sub1),]
Y <- as(t(Y), "sparseMatrix")
Matrix::writeMM(Y, file="../data/061923_D9_m2_Swiss_fibroblast_counts.mtx")


meta.sub2 <- meta |> filter(Mouse_ID == "080823_D9_m5_Swiss") |>
  filter(Tier1 %in% c("Fibroblast", "Epithelial"))


xy <- meta.sub2 |> select(x,y) |> as.matrix()
fit <- CurveFinderInteractive(xy)

save(fit, file="../data/080823_D9_m5_Swiss_fibroblast.RData")

# Save matrix

Y <- ad$raw[rownames(meta.sub2),]
Y <- as(t(Y), "sparseMatrix")
Matrix::writeMM(Y, file="../data/080823_D9_m5_Swiss_fibroblast_counts.mtx")



meta.sub3 <- meta |> filter(Mouse_ID == "080823_D9_m13_Swiss") |>
  filter(Tier1 %in% c("Fibroblast", "Epithelial"))



xy <- meta.sub3 |> select(x,y) |> as.matrix()
fit <- CurveFinderInteractive(xy)

save(fit, file="../data/080823_D9_m13_Swiss_fibroblast.RData")

# Save matrix

Y <- ad$raw[rownames(meta.sub3),]
Y <- as(t(Y), "sparseMatrix")
Matrix::writeMM(Y, file="../data/080823_D9_m13_Swiss_fibroblast_counts.mtx")


