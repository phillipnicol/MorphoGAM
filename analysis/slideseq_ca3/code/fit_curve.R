

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

l.o <- log(colSums(Y))
coef <- rep(0, nrow(Y))
p.val <- rep(1,nrow(Y))
peak <- coef
range <- coef
nz <- apply(Y, 1, function(x) sum(x != 0))
t <- fit$xyt$t
#f <- matrix(0,nrow=nrow(Y), ncol=length(t))
for(i in 1:nrow(Y)) {
  if(nz[i] < 10) {
    next
  }
  print(i)
  fit1 <- mgcv::gam(Y[i,] ~ s(t,bs="cr") + offset(l.o),family=nb(), H=diag(10))

  se_beta <- sqrt(diag(fit1$rV %*% t(fit1$rV))[-1])

  #Shrink
  my.ash <- ashr::ash(fit1$coefficients[-1], se_beta)
  beta.shrink <- apply(ashr::get_post_sample(my.ash,1000),
                       2,
                       median)

  mat <- as.matrix(mgcv::predict.gam(fit1, type = "lpmatrix")[,-1])
  fx1 <- as.vector(mat %*% beta.shrink)

  #fx.t[i,] <- fx1
  #fx.r[i,] <- fx2

  peak[i] <- max(fx1)
  p.val[i] <- summary(fit1)$s.pv
  range[i] <- max(exp(fit1$coefficients[1] + fx1)) - min(exp(fit1$coefficients[1] + fx1))
}

results_df <- data.frame(p.val = p.val, peak=peak, range=range)
rownames(results_df) <- rownames(Y)

saveRDS(results_df, file="../data/ca3_svg.RDS")

peak <- res$peak; range <- res$range
Y <- as.matrix(Y)
top5peak <- order(peak, decreasing=TRUE)[1:5]
expr <- t(Y[top5peak,]) |>
  as.data.frame() |>
  mutate(t=fit$xyt$t) |>
  pivot_longer(cols=-c(t)) |>
  mutate(name=fct_inorder(name))


ppeak <- ggplot(data=expr,aes(x=t,y=value)) +
  geom_point(size=0.5) +
  facet_wrap(~name, nrow=1, scales="free_y") +
  ylab("count") + theme_bw() +
  ggtitle("Peak")


top5range <- order(range, decreasing=TRUE)[1:5]
expr <- t(Y[top5range,]) |>
  as.data.frame() |>
  mutate(t=fit$xyt$t) |>
  pivot_longer(cols=-c(t)) |>
  mutate(name=fct_inorder(name))


prange <- ggplot(data=expr,aes(x=t,y=value)) +
  geom_point(size=0.5) +
  facet_wrap(~name, nrow=1, scales="free_y") +
  xlab("t") + ylab("count") + theme_bw() +
  ggtitle("Range")

library(ggpubr)
p <- ggarrange(ppeak, prange, nrow=2)



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
  geom_point(size=0.5, col="grey") +
  geom_line(aes(y=exp(fit1$coefficients[1] + fx1)), color="red") +
  scale_y_sqrt() +
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
  geom_point(size=0.5, col="grey") +
  geom_line(aes(y=exp(fit1$coefficients[1] + fx1.2)), color="red") +
  scale_y_sqrt() +
  theme_bw() +
  ylab("Count") +
  xlab("t") +
  ggtitle("Cpne9")

p.cable <- ggarrange(pgsr14, pcpne9,nrow=1)
ggsave(p.cable, filename="../plots/cable_genes.png",
       width=7.4, height=5.33, units="in")
