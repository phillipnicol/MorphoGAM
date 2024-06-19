

Y <- read.csv("../../data/GL2_distal_colon_cell_by_gene_raw.csv")

rownames(Y) <- Y[,1]

Y <- Y[,-1]

library(Matrix)

Y <- as.matrix(Y)

Y <- t(Y) #Transpose to get genes x cells


meta <- read.csv("../../data/GL2_distal_colon_cell_type_and_locations_2023.08.11.csv")

library(tidyverse)
meta <- as.data.frame(meta)
meta.sub <- meta |> filter(slice_full_name == "20220518_WT_dcol_slice_3") |>
  filter(spatial_neighborhood_v1 == "Mucosa") |>
  filter(leiden_combined_v2 == "Enterocyte")

Y.sub <- Y[,meta.sub$X]

Y.sub <- Y.sub[rowSums(Y.sub) >= 10,]

old.rownames <- rownames(Y.sub)

#Test expression for these genes
library(SPARK)
locus <- as.matrix(meta.sub[,c("x","y")])
res <- SPARK::sparkx(count_in = Y.sub,
                     locus_in = locus)

spark <- res$res_mtest |> as.data.frame() |>
  arrange(adjustedPval)

# nnSVG
library(nnSVG)
logCPM <- log(sweep(Y.sub, MARGIN=2, STATS = colSums(Y.sub), FUN="/")+1)
nn.svg <- nnSVG(input=logCPM, spatial_coords = locus,
                verbose=TRUE)

saveRDS(nn.svg, file="../data/nnSVG_results.RDS")

### Plot top 9


#Plot these two
expr <- t(Y.sub[rownames(nn.svg)[1:9],]) |>
  apply(2,function(x) log2(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))

nnsvgplot <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.25,alpha=0.75) +
  scale_color_gradient(low="grey90", high="darkred")+
  facet_wrap(~name) +
  labs(color="log expression") +
  theme_bw() + ggtitle("Top 9 SVGs (nnSVG)")
ggsave(nnsvgplot, filename="../plots/nnsvg_top9.png",width=13.4,
       height=10.1)

#Plot these two
expr <- t(Y.sub[rownames(spark)[1:9],]) |>
  apply(2,function(x) log2(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))
sparkplot <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.25,alpha=0.75) +
  scale_color_gradient(low="grey90", high="darkred")+
  facet_wrap(~name) +
  labs(color="log expression") +
  theme_bw() +
  ggtitle("Top 9 SVGs (SPARK)")
ggsave(sparkplot, filename="../plots/spark_top9.png",width=13.4,
       height=10.1)

library(mgcv)

radial.p <- rep(1, nrow(Y.sub))
angle.p <- rep(1, nrow(Y.sub))
angle <- rep(0,nrow(Y.sub))
radial <- rep(0,nrow(Y.sub))
angle.response <- angle
radial.response <- radial
l.o <- log(colSums(Y.sub))
fit <- CurveFinder(locus,knn=10,loop=TRUE)
fx.r <- Y.sub
fx.t <- Y.sub

for(i in 1:nrow(Y.sub)) {
  print(i)
  fit1 <- mgcv::gam(Y.sub[i,] ~ s(fit$xyt$t,bs="cc") + s(fit$xyt$r,bs="cr") + offset(l.o), family=nb(),
                    H = diag(18))

  se_beta <- sqrt(diag(fit1$rV %*% t(fit1$rV))[-1])

  #Shrink
  my.ash <- ashr::ash(fit1$coefficients[-1], se_beta)
  beta.shrink <- apply(ashr::get_post_sample(my.ash,1000),
                       2,
                       median)

  mat <- as.matrix(mgcv::predict.gam(fit1, type = "lpmatrix")[,-1])
  fx1 <- as.vector(mat[,1:8] %*% beta.shrink[1:8])
  fx2 <- as.vector(mat[,9:17] %*% beta.shrink[9:17])

  fx.t[i,] <- fx1
  fx.r[i,] <- fx2

  angle.response[i] <- max(exp(fit1$coefficients[1] + fx1) - exp(fit1$coefficients[1]))
  radial.response[i] <- max(exp(fit1$coefficients[1] + fx2) - exp(fit1$coefficients[1]))

  #pp <- predict(fit1, type="terms")

  angle[i] <- max(abs(fx1))
  angle.p[i] <- summary(fit1)$s.pv[1]
  radial[i] <- max(abs(fx2))
  radial.p[i] <- summary(fit1)$s.pv[2]
}

rownames(Y.sub) <- old.rownames
top5 <- order(angle, decreasing=TRUE)[1:5]
for(i in top5) {
  spark_rank <- which(rownames(p.vals) == rownames(Y.sub)[i])
  nnsvg_rank <- 0
  rownames(Y.sub)[i] <- paste("Spark rank =", spark_rank,
                              "nnSVG rank =", nnsvg_rank)
}

expr <- t(Y.sub[top5,]) |>
  as.data.frame() |>
  mutate(t=fit$xyt$t) |>
  pivot_longer(cols=-c(t))

p1 <- ggplot(data=expr,aes(x=t,y=value)) +
  geom_point() +
  facet_wrap(~name, ncol=2, scales="free_y")+
  xlab("t") + theme_bw()

top5 <- order(radial, decreasing=TRUE)[1:5]
for(i in top5) {
  spark_rank <- which(rownames(p.vals) == rownames(Y.sub)[i])
  nnsvg_rank <- 0
  rownames(Y.sub)[i] <- paste("Spark rank =", spark_rank,
                              "nnSVG rank =", nnsvg_rank)
}

expr <- t(Y.sub[top5,]) |>
  as.data.frame() |>
  mutate(t=fit$xyt$r) |>
  pivot_longer(cols=-c(t))

p2 <- ggplot(data=expr,aes(x=t,y=value)) +
  geom_point() +
  facet_wrap(~name, ncol=1, scales="free_y") +
  xlab("r") + theme_bw()


p.circ <- fit$curve.plot + guides(color="none")+
  xlab("") + ylab("") + ggtitle("Fitted curve")

p.resid <- fit$residuals.plot + guides(color="none") +
  xlab("") + ylab("")

ggsave(ggarrange(p.circ,p.resid,p1,p2,ncol=2,nrow=2,
                 heights=c(1,2)),
       filename="../plots/circle_genes.png",
       width=6, height=12)
