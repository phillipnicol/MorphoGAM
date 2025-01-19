library(tidyverse)
library(mgcv)

#results_df <- readRDS("../data/ca3_svg.RDS")

#peak <- results_df$peak; range <- results_df$range



library(STexampleData)

spe <- STexampleData::SlideSeqV2_mouseHPC()

ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]
Y <- counts(spe)[,ixs]

xy.dist <- as.matrix(dist(xy))
knn <- 20
prune.outlier <- 3
nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
outlier <- which(nnk > (prune.outlier)*median(nnk))

xy <- xy[-outlier,]; Y <- Y[,-outlier]

fit <- CurveFinder(xy)

load("../data/mgam_ca3.RData")


##First do the previous reported genes
Y <- as.matrix(Y)
p.prev <- plotGAMestimates(Y,
                           genes=c("Rgs14", "Cpne9"),
                           mgam_object = mgam,
                           curve_fit = fit) +
  scale_y_sqrt()





## Top 5 peak and range

p.peak <- plotGAMestimates(Y,
                           genes=order(mgam$results$peak.t,
                                       decreasing = TRUE)[1:5],
                           mgam_object = mgam,
                           curve_fit=fit) +
  scale_y_sqrt()

p.range <- plotGAMestimates(Y,
                           genes=order(mgam$results$range.t,
                                       decreasing = TRUE)[1:5],
                           mgam_object = mgam,
                           curve_fit=fit) +
  scale_y_sqrt()


## Sim power analysis
res <- readRDS("../data/results.RDS")
res <- res[,,-4]
library(reshape2)
library(tidyverse)
kappa <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
sigma <- seq(100, 1000, by=100)
df <- melt(res); colnames(df) <- c("kappa", "sigma", "method", "power")
df$kappa <- paste("k = ", kappa[df$kappa])
df$sigma <- sigma[df$sigma]

df$method <- c("SPARK-X", "MorphoGAM", "nnSVG")[df$method]

p.power <- ggplot(data=df,aes(x=sigma, y=power, color=method)) +
  geom_point() +
  geom_line() +
  facet_wrap(~kappa)+
  labs(color="Method")+
  xlab(expression(sigma)) + ylab("Power")+
  theme_bw() + theme(legend.position = "top")


p.curve <- fit$curve.plot + ggtitle("CA3 cells")
p.full <- ggarrange(ggarrange(p.curve, p.prev, nrow=1, ncol=2,
                              labels=c("a","b","")),
          p.power,
          ggarrange(p.peak,
                    p.range,
                    nrow=2),
          nrow=3,
          heights=c(1.5, 2, 2),
          labels=c("","c", "d"))
ggsave(p.full, filename="../plots/ca3_full_plot.png",
       width=9.11*1.1, height=11.7*1.1, units="in")



## For talk

p.curve.pow <- ggarrange(ggarrange(p.curve, p.prev, nrow=1, ncol=2),
                                   p.power, nrow=2,
                         heights=c(1.5,2))

ggsave(p.curve.pow, filename="../plots/ca3_curve_and_pow.png")

p.new.genes <- ggarrange(p.peak,
                         p.range,
                         nrow=2)

ggsave(p.new.genes, filename="../plots/ca3_new_genes.png")

## Supplement localization

df <- reshape2::melt(v) |>
  group_by(Var1, Var2, Var4) |>
  summarise(mean=mean(value))

p <- ggplot(data=df,aes(x=Var2, y=mean, color=Var4)) +
  geom_point() +
  facet_wrap(~Var1)
