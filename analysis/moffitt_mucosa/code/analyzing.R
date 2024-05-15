

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
logCPM <- log(sweep(Y.sub, MARGIN=2, STATS = colSums(Y.sub), FUN="/")+1)
nn.svg <- nnSVG(input=logCPM, spatial_coords = locus,
                verbose=TRUE)

library(mgcv)

radial.p <- rep(1, nrow(Y.sub))
angle.p <- rep(1, nrow(Y.sub))
angle <- rep(0,nrow(Y.sub))
radial <- rep(0,nrow(Y.sub))
l.o <- log(colSums(Y.sub))
fit <- CurveFinder(locus,knn=10,loop=TRUE)

fit1 <- mgcv::gam(Y.sub[1,] ~ s(fit$xyt$t,bs="cc") + offset(l.o), family=nb(),
                  fit=FALSE)
lambda1 <- createPenalty(fit1$X[,-1])
fit2 <- mgcv::gam(Y.sub[1,] ~ s(fit$xyt$r,bs="cr") + offset(l.o), family=nb(),
                  fit=FALSE)
lambda2 <- createPenalty(fit2$X[,-1])

for(i in 1:nrow(Y.sub)) {
  print(i)
  fit1 <- mgcv::gam(Y.sub[i,] ~ s(fit$xyt$t,bs="cc") + s(fit$xyt$r,bs="cr") + offset(l.o), family=nb(),
                    H = diag(18))
  pp <- predict(fit1, type="terms")
  angle[i] <- max(pp[,1])
  angle.p[i] <- summary(fit1)$s.pv[1]
  radial[i] <- max(abs(pp[,2]))
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
  facet_wrap(~name, ncol=1, scales="free_y")+
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
