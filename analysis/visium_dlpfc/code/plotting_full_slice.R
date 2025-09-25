


setwd(here::here("analysis", "visium_dlpfc", "code"))

library(tidyverse)
library(ggpubr)

load("../data/mgam_all.RData")

df <- mgam$results |> filter(pv.r < 0.05/(2*nrow(Y.input))) |> arrange(desc(peak.r)) |> head(n=12)

plotGAMestimates(Y.input, genes=rownames(df), curve_fit=fit, mgam_object = mgam, nrow=3,type="r")

df <- mgam$results |> arrange(desc(range.r)) |> head(n=10)


my.svd <- irlba::irlba(mgam$fxs.r[mgam$results$pv.r < 0.05/nrow(Y.input),], nv=10, nu=10)

Y.sub <- Y.input[mgam$results$pv.r < 0.05/nrow(Y.input),]

plotGAMestimates(Y.input, genes=c("SEMA3E", "CUX2"), curve_fit=fit, mgam_object = mgam, nrow=1,type="r")

plotGAMestimates(Y.input, genes=c("TRABD2A", "PENK"), curve_fit=fit, mgam_object = mgam, nrow=1,type="r")

plotGAMestimates(Y.input, genes=c("ADCYAP1", "ALB"), curve_fit=fit, mgam_object = mgam, nrow=1,type="r")


plotGAMestimates(Y.input, genes=rownames(df), curve_fit=fit, mgam_object = mgam, nrow=2,type="r")

plotFPCloading(mgam_object = mgam,curve.fit=fit,L=3,type="r")


plotGAMestimates(Y.input, genes=c("FABP7","PVALB", "CCK", "ENC1"), curve_fit=fit, mgam_object = mgam, nrow=1,type="r")

plotGAMestimates(Y.input, genes=c("AQP4","TRABD2A", "HPCAL1", "KRT17"), curve_fit=fit, mgam_object = mgam, nrow=1,type="r")

plotGAMestimates(Y.input, genes=c("AQP4","TRABD2A", "HPCAL1", "KRT17"), curve_fit=fit, mgam_object = mgam, nrow=1,type="r")


plotGAMestimates(Y.input, genes=c("MOBP","PCP4", "SNAP25", "KRT17"), curve_fit=fit, mgam_object = mgam, nrow=1,type="r")



df <- mgam$results |> filter(pv.r < 0.05/(2*nrow(Y.input))) |> arrange(desc(peak.r)) |> head(n=12)
df <- df |> arrange()
plotGAMestimates(Y.input, genes=rownames(df), curve_fit=fit, mgam_object = mgam, nrow=3,type="r")

df <- data.frame(x=fit$xyt$x,y=fit$xyt$y, color=my.svd$v[,3])
p <- ggplot(data=df,aes(x=x,y=y,color=color)) + geom_point()

gene <- "FOXP2"

ix <- which(rownames(counts(spe)) == "ENSG00000131095")
p.gene <- data.frame(x=spatialCoords(spe)[,1], y=spatialCoords(spe)[,2], color=counts(spe)[ix,]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient2(low="blue",mid="grey",high="red",midpoint=median(Y.input[gene,]))


plotGAMestimates(Y.input, genes=rownames(Y.sub)[order(my.svd$u[,3],decreasing=TRUE)[1:5]], curve_fit=fit, mgam_object = mgam, nrow=1,type="r")


p.existing_genes <- plotGAMestimates(Y.input, mgam_object=mgam,
                 curve_fit=fit,genes=c("AQP4", "TRABD2A", "HPCAL1", "KRT17"),nrow=2,type="r")


mgam$results[c("AQP4", "TRABD2A", "HPCAL1", "KRT17"),]
mgam$results["GFAP",]

gene.expr <- log(Y.input["GFAP",]/median(colSums(Y.input)) + 1)
p.gene <- data.frame(x=fit$xyt$x, y=fit$xyt$y, color=gene.expr) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient2(low="blue",mid="grey",high="red",midpoint=median(gene.expr)) +
  theme_bw() + ggtitle("GFAP log expression") + guides(color="none")

p.gene.fitted <- data.frame(x=fit$xyt$r, y=mgam$fxs.r["GFAP",]) |>
  ggplot(aes(x=x,y=y)) + geom_line() +
  theme_bw() + ggtitle("GFAP fitted function") + xlab("r") + ylab("Fitted function")


p.spatial <- ggarrange(p.layers +labs(color="ground truth"),p.t + guides(color="none") + ggtitle("t coordinate"),
                       p.r + guides(color="none") + ggtitle("r coordinate"),nrow=1,
                       labels=c("a","b",""))

p.all <- ggarrange(p.spatial, p.existing_genes,
                   ggarrange(p.gene.fitted,p.gene,nrow=1),
                   nrow=3, heights=c(1,1.3,1),
                   labels=c("","c","d"))

ggsave(p.all, filename="../plots/dlpca_all_slice.png",
       width=1.5*6.09, height=1.5*5.2, units="in")
