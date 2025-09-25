
setwd(here::here("analysis", "visium_dlpfc", "code"))

library(STexampleData)

library(MorphoGAM)

library(tidyverse)

library(biomaRt)

spe <- Visium_humanDLPFC()

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

Y.big <- counts(spe) |> as.matrix()

# Get mapping
genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters   = "ensembl_gene_id",
  values    = rownames(Y.big),
  mart      = ensembl
)

# keep only rows with non-empty symbols
genes <- genes[genes$hgnc_symbol != "" & !is.na(genes$hgnc_symbol), ]

# remove any symbols that occur more than once
dup_syms <- genes$hgnc_symbol[duplicated(genes$hgnc_symbol)]
genes <- genes[!genes$hgnc_symbol %in% dup_syms, ]

# make sure mapping is aligned to Y
genes <- genes[match(intersect(rownames(Y.big), genes$ensembl_gene_id),
                     genes$ensembl_gene_id), ]

# subset Y and rename
Y.big <- Y.big[genes$ensembl_gene_id, ]
rownames(Y.big) <- genes$hgnc_symbol




xy <- spatialCoords(spe)

ixs <- which(spe$ground_truth %in% c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6"))

p <- data.frame(x=xy[ixs,1],y=xy[ixs,2],color=spe$ground_truth[ixs]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point()

p.layers <- p

library(MorphoGAM)

xy <- xy[ixs,]

jxs <- which(spe$ground_truth[ixs] %in% c("Layer4"))

fit <- MorphoGAM::CurveFinder(xy[jxs,], k=10)

xyt <- fit$xyt |> arrange(t)
ft <- cbind(xyt$f1, xyt$f2)

proj <- princurve::project_to_curve(xy,s=as.matrix(ft))

t <- as.numeric(proj$lambda/max(proj$lambda))

t.full <- t

t <- fit$xyt$t

fitx <- mgcv::gam(xy[jxs,1]~s(t))
fity <- mgcv::gam(xy[jxs,2]~s(t))

all.fitx <- predict(fitx, newdata=data.frame(t=t.full))
all.fity <- predict(fity, newdata=data.frame(t=t.full))

rx <- xy[,1] - all.fitx
ry <- xy[,2] - all.fity


f2x <- gratia::derivatives(fitx,order=1,data=data.frame(t=t))
f2y <- gratia::derivatives(fity,order=1,data=data.frame(t=t))
t2 <- t

for(i in 1:length(t.full)) {
  print(i)
  e <- c(rx[i], ry[i])

  my.dist <- sqrt((all.fitx[i] - fitted(fitx))^2 + (all.fity[i] - fitted(fity))^2)
  t.anchor <- which.min(my.dist)

  Rf1 <- c(-f2y$.derivative[t.anchor], f2x$.derivative[t.anchor])
  sign <- ifelse(sum(e*Rf1) > 0, 1, -1)
  t2[i] <- sign*sqrt(sum(e^2))
}

t2 <- (t2 - min(t2))/(max(t2) - min(t2))

r <- t2


## Smooth t in one final step

smooth.t <- mgcv::gam(t.full ~ s(xy[,1], xy[,2]))



p <- data.frame(x=xy[,1],y=xy[,2],color=fitted(smooth.t)) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient2(low="blue",mid="grey",high="red",midpoint=0.5)

p.t <- p

p <- data.frame(x=xy[,1],y=xy[,2],color=r) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient2(low="blue",mid="grey",high="red",midpoint=0.5)

p.r <- p

library(ggplot2)
library(ggpubr)
p <- ggarrange(p.layers,p.t,p.r,nrow=1)

ggsave(p, filename="../plots/dlpfc_layers_morpho.png", width=9, height=3)


fit <- list()
fit$xyt <- data.frame(x=xy[,1],y=xy[,2],t=t.full, r=r)

Y.input <- Y.big[,ixs]
Y.input <- Y.input[rowSums(Y.input) >= 25,]

#mgam <- MorphoGAM(Y=Y.input[1:100,],curve.fit=fit,
#                  design=y ~ s(t, bs="cr", k=10) + s(r, bs="cr", k=10))

mgam <- MorphoGAM(Y=Y.input,curve.fit=fit,
                  design = y ~ s(t, bs="cr", k=10) + s(r, bs ="cr", k=10))


save(mgam, file="../data/mgam_all.RData")

print("COMPLETED")
