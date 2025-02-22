setwd(here::here("analysis/merfish_swissroll/code"))

library(tidyverse)
library(ggpubr)

#library(MorphoGAM)

#xy <- matrix(cbind(ad$obs$x, ad$obs$y),ncol=2)


#fit <- CurveFinderInteractive(xy)
#
#save(fit, file="../data/curve_d9_061923_baysor.RData")

load("../data/mgam_d9_061923_baysor.RData")

load("../data/curve_d9_061923.RData")

gene_names <- readRDS("../data/gene_names.RDS")

library(MorphoGAM)

Y <- Matrix::readMM("../data/d9_061923.mtx")
Y <- Matrix::t(Y)

rownames(Y) <- gene_names
rownames(mgam$results) <- gene_names
rownames(mgam$fxs.r) <- gene_names
rownames(mgam$fxs.t) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names


MorphoGAM::plotGAMestimates(Y,genes=c("Mstn", "Clu", "Ltb", "Cd22"),
                            type="r", curve_fit=fit, mgam_object=mgam, nrow=2)

mgam1 <- mgam
fit1 <- fit


load("../data/curve_d9_m5_080823.RData")
Y <- Matrix::readMM("../data/080823_d9_m5.mtx")
load("../data/mgam_d9_m5_080823_baysor.RData")

Y <- Matrix::t(Y)

rownames(Y) <- gene_names
rownames(mgam$results) <- gene_names
rownames(mgam$fxs.r) <- gene_names
rownames(mgam$fxs.t) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names

mgam2 <- mgam
fit2 <- fit





load("../data/curve_d9_m13_080823.RData")
Y <- Matrix::readMM("../data/080823_d9_m13.mtx")
load("../data/mgam_d9_m13_080823_baysor.RData")

Y <- Matrix::t(Y)

rownames(Y) <- gene_names
rownames(mgam$results) <- gene_names
rownames(mgam$fxs.r) <- gene_names
rownames(mgam$fxs.t) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names

mgam3 <- mgam
fit3 <- fit


meta <- data.frame(peak.t.1 = mgam1$results$peak.t,
                   range.t.1 = mgam1$results$range.t,
                   peak.r.1= mgam1$results$peak.r,
                   range.r.1=mgam1$results$range.r,
                   peak.t.2 = mgam2$results$peak.t,
                   range.t.2 = mgam2$results$range.t,
                   peak.r.2= mgam2$results$peak.r,
                   range.r.2=mgam2$results$range.r,
                   peak.t.3 = mgam3$results$peak.t,
                   range.t.3 = mgam3$results$range.t,
                   peak.r.3= mgam3$results$peak.r,
                   range.r.3=mgam3$results$range.r)

rownames(meta) <- gene_names

peak.t.min <- apply(meta,1,function(x) min(x[c(1,4,7)]))

peak.t.min <- peak.t.min[order(peak.t.min, decreasing = TRUE)[1:6]]


reg.coef1 <- apply(mgam1$fxs.t, 1, function(x) {
  my.fit <- lm(x ~ fit1$xyt$t)
  my.fit$coefficient[2]
})

reg.coef2 <- apply(mgam2$fxs.t, 1, function(x) {
  my.fit <- lm(x ~ fit2$xyt$t)
  my.fit$coefficient[2]
})

reg.coef3 <- apply(mgam3$fxs.t, 1, function(x) {
  my.fit <- lm(x ~ fit3$xyt$t)
  my.fit$coefficient[2]
})

max.min.reg <- pmax(reg.coef1, reg.coef2, reg.coef3)

min.max.reg <- pmin(reg.coef1, reg.coef2, reg.coef3)


plotting.genes <- c("Spdef", "Cldn15", "Muc4", "Igfbp5", "C3", "Il21r")
df <- data.frame(t=c(), curve=c(), y=c(), gene=c())

for(gene in plotting.genes) {
  df <- rbind(df,data.frame(t=fit1$xyt$t, curve="Swissroll 1", y=mgam1$fxs.t[gene,],
                            gene=gene),
              data.frame(t=fit2$xyt$t, curve="Swissroll 2", y=mgam2$fxs.t[gene,],
                         gene=gene),
              data.frame(t=fit3$xyt$t, curve="Swissroll 3", y=mgam3$fxs.t[gene,],gene=gene))

}


p.genes <- ggplot(data=df,aes(x=t,y=y,color=curve)) + geom_line() +
  facet_wrap(~gene, nrow=2,ncol=3, scales="free_y") +
  geom_abline(slope=0, intercept=0, color="grey", linetype="dashed") +
  scale_color_manual(labels=c("Swissroll 1", "Swissroll 2", "Swissroll 3"),
                     values=c("darkred", "darkblue", "darkgreen")) +
  theme_bw() + labs(color="") + ylab("Log FC from baseline") +
  theme(legend.position = "bottom")

p.swiss <- ggarrange(fit1$curve.plot + ggtitle("Swissroll 1") + labs(color="t"),
          fit2$curve.plot + ggtitle("Swissroll 2") + labs(color="t"),
          fit3$curve.plot + ggtitle("Swissroll 3") + labs(color="t"),
          nrow=1, common.legend = TRUE, legend="bottom")

p <- ggarrange(p.swiss, p.genes, nrow=2, labels=c("a","b"),
               heights=c(1,1.25))

ggsave(p, filename="../plots/swissroll_all_genes.png",
       width=9.44, height=7.78)


range.t.min <- apply(meta,1,function(x) min(x[c(2,5,8)]))
range.t.min <- range.t.min[order(range.t.min, decreasing = TRUE)[1:6]]



#Plot

df <- data.frame(x=fit2$xyt$x,y=fit2$xyt$y,
                 color=Y["Cldn15",]) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.25) +
  scale_color_gradient(low="grey90", high="red") +
  theme_bw()



