

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

Y.s <- sweep(Y.sub, MARGIN = 2, STATS = colSums(Y.sub)/median(colSums(Y.sub)),
             FUN="/")

expr <- t(Y.s[rownames(nn.svg)[1:9],]) |>
  apply(2,function(x) log2(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))

nnsvgplot <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.25) +
  scale_color_gradient(low="grey85", high="darkred")+
  #scale_size_continuous(range=c(0.1,0.75)) +
  scale_alpha_continuous(range=c(0.33,1)) +
  coord_fixed() +
  facet_wrap(~name) +
  labs(color="log expression") +
  theme_bw() +
  ggtitle("Top 9 SVGs (nnSVG)")
ggsave(nnsvgplot, filename="../plots/nnsvg_top9.png",
       width=1.3*6.42, height=1.3*5.85,
       units="in")

#Plot these two
expr <- t(Y.s[rownames(spark)[1:9],]) |>
  apply(2,function(x) log2(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))
sparkplot <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.25) +
  scale_color_gradient(low="grey85", high="darkred")+
  #scale_size_continuous(range=c(0.1,0.75)) +
  scale_alpha_continuous(range=c(0.33,1)) +
  coord_fixed() +
  facet_wrap(~name) +
  labs(color="log expression") +
  theme_bw() +
  ggtitle("Top 9 SVGs (SPARKX)")
ggsave(sparkplot, filename="../plots/spark_top9.png",
       width=1.3*6.42, height=1.3*5.85,
       units="in")

## MorphoGAM
library(mgcv)
fit <- CurveFinder(locus,knn=10,loop=TRUE)

mgam <- MorphoGAM(Y.sub,
                  curve.fit=fit,
                  design=y~s(t,bs="cc")+s(r,bs="cr"))


save(mgam, file="../data/mucosa_mgam.RData")
save(fit, file="../data/curve.RData")


#results_df <- readRDS("../data/loop_svg.RDS")
nnsvg_res <- readRDS("../data/nnSVG_results.RDS")

#angle <- results_df$angle
#radial.response <- results_df$radial.response
#nnsvg_rank <- nnsvg_res$rank
#p.vals <- spark$adjustedPval

load("../data/mucosa_mgam.RData")
load("../data/curve.RData")

rownames(Y.sub) <- old.rownames
top5 <- order(mgam$results$peak.t, decreasing=TRUE)[1:6]
#top5 <- result[1:5]
for(i in top5) {
  spark_rank <- which(rownames(spark) == rownames(Y.sub)[i])
  nnsvg_rank <- which(rownames(nnsvg_res) == rownames(Y.sub)[i])
  rownames(Y.sub)[i] <- paste0(old.rownames[i], ": ",
                               "Spark rank = ", spark_rank,
                               " nnSVG rank = ", nnsvg_rank)
}


p.peak <- plotGAMestimates(Y.sub,
                           genes=order(mgam$results$peak.t,decreasing = TRUE)[1:6],
                           curve_fit = fit,
                           mgam_object = mgam,
                           nrow=2) +
  theme(strip.text=element_text(size=3.5*1.5)) +
  scale_y_sqrt()

rownames(Y.sub) <- old.rownames
top5 <- order(mgam$results$range.r, decreasing=TRUE)[1:6]
#top5 <- result[1:5]
for(i in top5) {
  spark_rank <- which(rownames(spark) == rownames(Y.sub)[i])
  nnsvg_rank <- which(rownames(nnsvg_res) == rownames(Y.sub)[i])
  rownames(Y.sub)[i] <- paste0(old.rownames[i], ": ",
                               "Spark rank = ", spark_rank,
                               " nnSVG rank = ", nnsvg_rank)
}

p.range <- plotGAMestimates(Y.sub,
                           genes=order(mgam$results$range.r,decreasing = TRUE)[1:6],
                           curve_fit = fit,
                           mgam_object = mgam,
                           type="r",
                           nrow=2) +
  theme(strip.text=element_text(size=3.5*1.5)) +
  xlim(0.3,0.7) +
  scale_y_sqrt()

library(ggpubr)
p <- ggarrange(p.peak, p.range, nrow=2, labels=c("a","b"))

ggsave(p, filename="../plots/mouse_mucosa_svgs.png",
       width= 4.78*1.5,
       height= 5.18*1.5)


mgam$results["Ddx58",]

mgam$results["Apob",]




rownames(Y.sub) <- old.rownames
top5 <- order(mgam$results$range.t, decreasing=TRUE)[1:6]
#top5 <- result[1:5]
for(i in top5) {
  spark_rank <- which(rownames(spark) == rownames(Y.sub)[i])
  nnsvg_rank <- which(rownames(nnsvg_res) == rownames(Y.sub)[i])
  rownames(Y.sub)[i] <- paste0(old.rownames[i], ": ",
                               "Spark rank = ", spark_rank,
                               " nnSVG rank = ", nnsvg_rank)
}


p.peak <- plotGAMestimates(Y.sub,
                           genes=order(mgam$results$range.t,decreasing = TRUE)[1:6],
                           curve_fit = fit,
                           mgam_object = mgam,
                           nrow=2) +
  theme(strip.text=element_text(size=3.5*1.5)) +
  scale_y_sqrt()

rownames(Y.sub) <- old.rownames
top5 <- order(mgam$results$peak.r, decreasing=TRUE)[1:6]
#top5 <- result[1:5]
for(i in top5) {
  spark_rank <- which(rownames(spark) == rownames(Y.sub)[i])
  nnsvg_rank <- which(rownames(nnsvg_res) == rownames(Y.sub)[i])
  rownames(Y.sub)[i] <- paste0(old.rownames[i], ": ",
                               "Spark rank = ", spark_rank,
                               " nnSVG rank = ", nnsvg_rank)
}

p.range <- plotGAMestimates(Y.sub,
                            genes=order(mgam$results$peak.r,decreasing = TRUE)[1:6],
                            curve_fit = fit,
                            mgam_object = mgam,
                            type="r",
                            nrow=2) +
  theme(strip.text=element_text(size=3.5*1.5)) +
  xlim(0.3,0.7) +
  scale_y_sqrt()

library(ggpubr)
p <- ggarrange(p.peak, p.range, nrow=2, labels=c("a","b"))

ggsave(p, filename="../plots/mouse_mucosa_svgs_OTHER.png",
       width= 4.78*1.5,
       height= 5.18*1.5)




