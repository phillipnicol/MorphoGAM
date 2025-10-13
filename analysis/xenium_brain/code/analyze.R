
setwd(here::here("analysis", "xenium_brain", "code"))


Y <- Matrix::readMM("../data/matrix.mtx")

genes <- read.csv("../data/features.tsv", sep="\t", header=FALSE)

my.cells <- read.csv("../data/cells.csv")

ix <- sample(1:nrow(my.cells), size=10^4)

meta.1 <- my.cells[ix,]

library(readxl)
ct <- read_excel("../data/41467_2023_44560_MOESM4_ESM.xlsx",
                                           sheet = "Fig5a")

library(tidyverse)


ixs <- which(ct$scClassify %in% c("CA1", "CA2", "CA3"))
cx <- rep(0, length(ixs))

for(i in 1:length(ixs)) {
  print(i)
  my.dist <- sqrt((my.cells$x_centroid - ct$cell_centroid_x[ixs[i]])^2 + (my.cells$y_centroid - ct$cell_centroid_y[ixs[i]])^2)
  cx[i] <- which.min(my.dist)
}


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

Y <- as.matrix(Y)

rownames(Y) <- genes$V2

t.old <- fit$xyt$t

#fit$xyt$t <- 2*abs(fit$xyt$t - 0.5) #Transform

fit$xyt$t.old <- t.old

mgam <- MorphoGAM(Y, curve.fit=fit,
                  design=y~s(t,bs="cr") + s(r,bs="cr"))


rownames(mgam$results) <- genes$V2

mgam$results |> arrange(desc(peak.t)) |> head()

plotGAMestimates(Y, genes=c("Strip2", "Syndig1", "Prox1", "Cpne8"), curve_fit=fit, mgam_object = mgam)

mgam$results |> arrange(desc(range.r)) |> head()

plotGAMestimates(Y, genes=c("Cabp7", "Epha4", "Neurod6", "Nrn1"), curve_fit=fit, mgam_object = mgam, type="r")


plotGAMestimates(Y, genes=c("Gfap", "Bcl11b", "Fibcd1", "Igfbp4"), curve_fit=fit, mgam_object = mgam)

expr <- Y["Pou3f1",]

#plot expression of gfap
df.expr <- data.frame(x=xy[,1], y=xy[,2], expr=expr)

p.expr <- df.expr |> ggplot(aes(x=x,y=y,color=expr)) + geom_point() +
  scale_color_viridis_c() + theme_bw() + labs(color="Gfap expression")




my.genes.r <- mgam$results |> arrange(desc(range.r)) |> head() |> rownames()

plotGAMestimates(Y, genes=my.genes.r, curve_fit=fit, mgam_object = mgam, type="r")


#plot expression of gfap
expr <- Y["Dkk3",]
df.expr <- data.frame(x=xy[,1], y=xy[,2], expr=(log(expr) + 1))

p.expr <- df.expr |> ggplot(aes(x=x,y=y,color=expr)) + geom_point() +
  scale_color_viridis_c() + theme_bw() + labs(color="Gfap expression")





#Test expression for these genes
library(SPARK)
locus <- as.matrix(xy)
res <- SPARK::sparkx(count_in = Y,
                     locus_in = locus)

spark <- res$res_mtest |> as.data.frame() |>
  arrange(adjustedPval)

saveRDS(spark, file="../data/spark_results.RDS")

plotGAMestimates(Y, genes=c("Pou3f1","Fibcd1","Ndst4","Satb2","Wfs1"), curve_fit=fit, mgam_object = mgam)
plotGAMestimates(Y, genes=c("Cpne4","Slit2","Cdh9","Rspo2","Galnt14"), curve_fit=fit, mgam_object = mgam)

fpca.t <- irlba::irlba(mgam$fxs.t[1:248,], nv=20, nu=20)

expr <- Y["Syt17",]
df.expr <- data.frame(x=xy[,1], y=xy[,2], expr=(log(expr) + 1))

p.expr <- df.expr |> ggplot(aes(x=x,y=y,color=expr)) + geom_point() +
  theme_bw() + labs(color="Gfap expression")


##Plots

fit$curve.plot

fit$coordinate.plot

fit$residuals.plot

p.firstpc <- plotGAMestimates(Y, genes=c("Pou3f1","Fibcd1","Ndst4","Satb2","Wfs1"), curve_fit=fit, mgam_object = mgam)
p.secondpc <- plotGAMestimates(Y, genes=c("Cpne4","Slit2","Cdh9","Rspo2","Galnt14"), curve_fit=fit, mgam_object = mgam)


write.csv(xy,"../../data/xenium_brain_coords.csv")
write.csv(Y[1:248,],"../../data/xenium_brain_counts.csv")
write.csv(rownames(Y)[1:248],"../../data/xenium_brain_gene_names.csv")


### Plotting


xyt <- fit$xyt |> arrange(desc(t.old))
df.new <- data.frame(x=xyt$x,
                     y=xyt$y)

p <- ggplot(data=df.new,aes(x=x,y=y)) + geom_point(col="grey",
                                                   alpha=0.5)

df.line <- data.frame(x=xyt$f1, y=xyt$f2)
p <- p + geom_path(data=df.line,aes(x=x,y=y),
                   linewidth=1)
p <- p + theme_bw() + ggtitle("CA Cells")


p2 <- data.frame(x=xyt$x,y=xyt$y,color=xyt$t) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),
                        colors=c("navyblue","grey90", "firebrick1"))+
  theme_bw() +
  ggtitle("First Coordinate") +
  labs(color="t") + guides(color="none")

p3 <- data.frame(x=xyt$x,y=xyt$y,color=xyt$r) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),
                        colors=c("navyblue","grey90", "firebrick1"))+
  theme_bw() +
  ggtitle("Second Coordinate") +
  labs(color="r") + guides(color="none")


library(ggpubr)

p.curves <- ggarrange(p, p2, p3, nrow=1)

p.firstpc  <- p.firstpc + ylab("")
p.secondpc <- p.secondpc + ylab("")
p.genes <- ggarrange(p.firstpc, p.secondpc, nrow=2, ncol=1)


p <- ggarrange(p.curves, p.genes, nrow=2, labels=c("a","b"), heights=c(1,1.5))

ggsave(p,filename="../plots/xenium_brain_fig.png",
       width=1.5*5.5, height=1.5*3.32, units="in")

