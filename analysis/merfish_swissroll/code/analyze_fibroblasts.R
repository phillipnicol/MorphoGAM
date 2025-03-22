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

# Find genes correlated with Vil1
vil1.1 <- mgam$fxs.t["Vil1",]

cos.sim1 <- apply(mgam$fxs.t, 1, function(x) {
  sum(x*vil1.1)/(sqrt(sum(x^2)*sum(vil1.1^2)))
})

mgam1 <- mgam


## Two


load("../data/mgam_d9_m5_080823_fibroblast_baysor.RData")

load("../data/080823_D9_m5_Swiss_fibroblast.RData")

Y <- Matrix::readMM("../data/080823_D9_m5_Swiss_fibroblast_counts.mtx")


rownames(Y) <- gene_names
Y2 <- Y
rownames(mgam$results) <- gene_names
rownames(mgam$fxs.r) <- gene_names
rownames(mgam$fxs.t) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names

# Find genes correlated with Vil1
vil1.2 <- mgam$fxs.t["Vil1",]

cos.sim2 <- apply(mgam$fxs.t, 1, function(x) {
  sum(x*vil1.2)/(sqrt(sum(x^2)*sum(vil1.2^2)))
})

cos.sim2[order(abs(cos.sim2), decreasing = TRUE)[1:10]]

mgam2 <- mgam


## Third



load("../data/mgam_d9_m13_080823_fibroblast_baysor.RData")

load("../data/080823_D9_m13_Swiss_fibroblast.RData")

Y <- Matrix::readMM("../data/080823_D9_m13_Swiss_fibroblast_counts.mtx")




rownames(Y) <- gene_names
Y3 <- Y
rownames(mgam$results) <- gene_names
rownames(mgam$fxs.r) <- gene_names
rownames(mgam$fxs.t) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names

# Find genes correlated with Vil1
vil1.3 <- mgam$fxs.t["Vil1",]

cos.sim3 <- apply(mgam$fxs.t, 1, function(x) {
  sum(x*vil1.3)/(sqrt(sum(x^2)*sum(vil1.3^2)))
})

cos.sim3[order(abs(cos.sim3), decreasing = TRUE)[1:10]]

mgam3 <- mgam


df <- data.frame(roll1 = cos.sim1,
                 roll2 = cos.sim2,
                 roll3 = cos.sim3)


df$neg.markers <- pmin(df$roll1, df$roll2, df$roll3)
df$pos.markers <- pmax(df$roll1, df$roll2, df$roll3)

df |> arrange(desc(neg.markers)) |>
  select(roll1, roll2, roll3) |> head(n=10)

df |> arrange(pos.markers) |>
  select(roll1, roll2, roll3) |> head(n=10)






