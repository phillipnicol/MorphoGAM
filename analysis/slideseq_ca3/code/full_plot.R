library(tidyverse)
library(mgcv)

results_df <- readRDS("../data/ca3_svg.RDS")

peak <- results_df$peak; range <- results_df$range



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


load("../data/ca3_curve.RData")

Y <- as.matrix(Y)

top5peak <- order(peak, decreasing=TRUE)[1:5]
expr <- t(Y[top5peak,]) |>
  as.data.frame() |>
  mutate(t=fit$xyt$t) |>
  pivot_longer(cols=-c(t)) |>
  mutate(name=fct_inorder(name))


ppeak <- ggplot(data=expr,aes(x=t,y=value)) +
  geom_jitter(width=0,height=0.1,size=0.5) +
  facet_wrap(~name, nrow=1, scales="free_y") +
  ylab("count") + theme_bw() +
  scale_y_continuous(trans="log1p") +
  ggtitle("Peak")


top5range <- order(range, decreasing=TRUE)[1:5]
expr <- t(Y[top5range,]) |>
  as.data.frame() |>
  mutate(t=fit$xyt$t) |>
  pivot_longer(cols=-c(t)) |>
  mutate(name=fct_inorder(name))


prange <- ggplot(data=expr,aes(x=t,y=value)) +
  geom_jitter(width=0,height=0.1,size=0.5) +
  facet_wrap(~name, nrow=1, scales="free_y") +
  xlab("t") + ylab("count") + theme_bw() +
  ggtitle("Range")

library(ggpubr)
p <- ggarrange(ppeak, prange, nrow=2)


l.o <- log(colSums(Y))
t <- fit$xyt$t
## Rgs14
i <- which(rownames(Y) == "Rgs14")
fit1 <- mgcv::gam(Y[i,] ~ s(t,bs="cr") + offset(l.o),family=nb(), H=diag(10))
print(fit1$sig2)

se_beta <- sqrt(diag(fit1$rV %*% t(fit1$rV))[-1])

#Shrink
my.ash <- ashr::ash(fit1$coefficients[-1], se_beta)
beta.shrink <- apply(ashr::get_post_sample(my.ash,1000),
                     2,
                     median)

mat <- as.matrix(mgcv::predict.gam(fit1, type = "lpmatrix")[,-1])
fx1 <- as.vector(mat %*% beta.shrink)

pgsr14 <- data.frame(x=fit$xyt$t, y=Y[i,]/exp(l.o)) |>
  ggplot(aes(x=x,y=y)) +
  geom_jitter(size=0.5, col="grey", width=0, height=0.00005) +
  geom_line(aes(y=exp(fit1$coefficients[1] + fx1)), color="red") +
  scale_y_continuous(trans="log1p") +
  theme_bw() +
  ylab("Count") +
  xlab("t") +
  ggtitle("Rgs14")

#Cpne9
i <- which(rownames(Y) == "Cpne9"); print(i)
fit1 <- mgcv::gam(Y[i,] ~ s(t,bs="cr") + offset(l.o),family=nb(), H=diag(10))
print(fit1$sig2)

se_beta <- sqrt(diag(fit1$rV %*% t(fit1$rV))[-1])

#Shrink
my.ash <- ashr::ash(fit1$coefficients[-1], se_beta)
beta.shrink <- apply(ashr::get_post_sample(my.ash,1000),
                     2,
                     median)

mat <- as.matrix(mgcv::predict.gam(fit1, type = "lpmatrix")[,-1])
fx1.2 <- as.vector(mat %*% beta.shrink)

pcpne9 <- data.frame(x=fit$xyt$t, y=Y[i,]/exp(l.o)) |>
  ggplot(aes(x=x,y=y)) +
  geom_jitter(size=0.5, col="grey", width=0, height=0.00005) +
  geom_line(aes(y=exp(fit1$coefficients[1] + fx1.2)), color="red") +
  scale_y_continuous(trans="log1p") +
  theme_bw() +
  ylab("Count") +
  xlab("t") +
  ggtitle("Cpne9")

p.cable <- ggarrange(pgsr14, pcpne9,nrow=1)

p.curve <- fit$curve.plot



res <- readRDS("../data/results.RDS")
res <- res[,,-4]
library(reshape2)
library(tidyverse)
kappa <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
sigma <- seq(100, 1000, by=100)
df <- melt(res); colnames(df) <- c("kappa", "sigma", "method", "power")
df$kappa <- paste("k = ", kappa[df$kappa])
df$sigma <- sigma[df$sigma]

df$method <- c("SPARK-X", "Projection", "nnSVG")[df$method]

p.power <- ggplot(data=df,aes(x=sigma, y=power, color=method)) +
  geom_point() +
  geom_line() +
  facet_wrap(~kappa)+
  labs(color="Method")+
  xlab(expression(sigma)) + ylab("Power")+
  theme_bw() + theme(legend.position = "top")


p.full <- ggarrange(ggarrange(p.curve, pgsr14, pcpne9, nrow=1, ncol=3,
                              labels=c("a","b","")),
          p,
          p.power,
          nrow=3,
          heights=c(1.5, 2, 2),
          labels=c("","d", "e"))
ggsave(p.full, filename="../plots/ca3_full_plot.png",
       width=9.11, height=11.7, units="in")
