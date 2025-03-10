setwd(here::here("analysis/merfish_swissroll/code"))


load("../data/mgam_d9_061923_fibroblast_baysor.RData")

load("../data/061923_D9_m2_Swiss_fibroblast.RData")

Y <- Matrix::readMM("../data/061923_D9_m2_Swiss_fibroblast_counts.mtx")

gene_names <- readRDS("../data/gene_names.RDS")

library(MorphoGAM)


rownames(Y) <- gene_names
Y1 <- Y
rownames(mgam$results) <- gene_names
rownames(mgam$fxs.r) <- gene_names
rownames(mgam$fxs.t) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names

t.range <- c(0.95, 1)
ixs <- which(fit$xyt$t > t.range[1] & fit$xyt$t < t.range[2])

my.mean <- rowMeans(mgam$fxs.t[,ixs])




