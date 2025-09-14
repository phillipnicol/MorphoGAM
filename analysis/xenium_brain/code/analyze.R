
setwd(here::here("analysis", "xenium_brain", "code"))


Y <- Matrix::readMM("../data/matrix.mtx")

genes <- read.csv("../data/features.tsv", sep="\t", header=FALSE)

my.cells <- read.csv("../data/cells.csv")

ix <- sample(1:nrow(my.cells), size=10^4)

meta.1 <- my.cells[ix,]

library(readxl)
ct <- read_excel("../data/41467_2023_44560_MOESM4_ESM.xlsx",
                                           sheet = "Fig5a")

ca_id <- vapply(ct$...1, FUN.VALUE=character(1), function(x) {
  str_split(x, "_")[[1]][2]
}) |> as.numeric()

library(tidyverse)

bidcell.sum <- ct |> group_by(cell_id) |> summarize(x=mean(coord_x),
                                                    y=mean(coord_y))

#ca_cells <- bidcell.sum[bidcell.sum$cell_type %in% c("CA1", "CA2", "CA3"),]

ca_cells <- bidcell.sum

meta.sub <- merge(my.cells, ca_cells, by=c("cell_id"="cell_id"))

p <- meta.sub |> mutate(x=x_centroid, y=y_centroid, color=cell_type) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point()



p <- meta.1 |> mutate(x=x_centroid, y=y_centroid) |>
  ggplot(aes(x=x,y=y)) + geom_point()

library(MorphoGAM)



1061.8757	1074.2216

ixs <- which(ct$scClassify %in% c("CA1", "CA2", "CA3"))
cx <- rep(0, length(ixs))

for(i in 1:length(ixs)) {
  print(i)
  cx[i] <- which(abs(my.cells$x_centroid - ct$cell_centroid_x[ixs[i]]) < 10 & abs(my.cells$y_centroid - ct$cell_centroid_y[ixs[i]]) < 10)[1]
}

cx <- cx[!is.na(cx)]

cx <- cx[which(my.cells$y_centroid[cx] > 4000)]

xy <- cbind(my.cells$x_centroid[cx], my.cells$y_centroid[cx])
xy.dist <- as.matrix(dist(xy))
knn <- 10
prune.outlier <- 3
nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
outlier <- which(nnk > (prune.outlier)*median(nnk))

xy <- xy[-outlier,]

library(MorphoGAM)

fit <- CurveFinder(xy,knn=10)

Y <- Y[,cx]
Y <- Y[,-outlier]

#Y <- t(as.matrix(Y))

mgam <- MorphoGAM(Y, curve.fit=fit,
                  design=y~s(t,bs="cr") + s(r,bs="cr"))

rownames(mgam$results) <- genes$V2
rownames(mgam$fxs.r) <- genes$V2
mgam$results |> arrange(desc(peak.t)) |> head()

plotGAMestimates(Y, genes=c("Cpne4", "Pou3f1", "Meis2", "Fibcd1"), curve_fit=fit, mgam_object = mgam)

plotGAMestimates(Y, genes=c("Gfap", "Pdgfra", "Unc13c", "Cldn5"), curve_fit=fit, mgam_object = mgam, type="r")


expr <- Y["Pou3f1",]

#plot expression of gfap
df.expr <- data.frame(x=xy[,1], y=xy[,2], expr=expr)

p.expr <- df.expr |> ggplot(aes(x=x,y=y,color=expr)) + geom_point() +
  scale_color_viridis_c() + theme_bw() + labs(color="Gfap expression")

