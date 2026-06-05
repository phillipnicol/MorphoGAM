setwd(here::here("analysis/moffitt_mucosa/code"))

library(RColorBrewer)
library(ggrepel)
library(ggpubr)

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

locus <- as.matrix(meta.sub[,c("x","y")])



## Some different choices of k plus hand drawn 
ktry <- c(10, 20, 50)

fit10 <- CurveFinder(locus,knn=10,loop=TRUE)
mgam10 <- MorphoGAM(Y.sub,
                  curve.fit=fit10,
                  design=y~s(t,bs="cc")+s(r,bs="cr"),
                  shrinkage = FALSE)


fit20 <- CurveFinder(locus,knn=20,loop=TRUE)
mgam20 <- MorphoGAM(Y.sub,
                  curve.fit=fit20,
                  design=y~s(t,bs="cc")+s(r,bs="cr"),
                  shrinkage = FALSE)


fit50 <- CurveFinder(locus,knn=50,loop=TRUE)
mgam50 <- MorphoGAM(Y.sub,
                  curve.fit=fit50,
                  design=y~s(t,bs="cc")+s(r,bs="cr"),
                  shrinkage = FALSE)

#Hand drawn 
load("../../principal_curves/data/loop_ground_truth.Rda")
mgam_hd <- MorphoGAM(Y.sub,
                  curve.fit=fit,
                  design=y~s(t,bs="cc")+s(r,bs="cr"),
                  shrinkage = FALSE)

#Save all the mgam
save(mgam10, file="../data/mucosa_mgam10.RData")
save(mgam20, file="../data/mucosa_mgam20.RData")
save(mgam50, file="../data/mucosa_mgam50.RData")
save(mgam_hd, file="../data/mucosa_mgam_hd.RData")

#Load them back
load("../data/mucosa_mgam10.RData")
load("../data/mucosa_mgam20.RData")
load("../data/mucosa_mgam50.RData")

#Load the one used in manuscript
load("../data/mucosa_mgam.RData")

#Plotting

library(ggplot2)
library(patchwork)

# Load the hand-drawn curve used in the manuscript
load("../data/mucosa_mgam.RData")

make_df <- function(ref, comp, stat, method) {
  data.frame(
    reference = ref,
    comparison = comp,
    stat = stat,
    method = method
  )
}

df <- rbind(
  make_df(mgam$results$peak.t, mgam10$results$peak.t, "t", "k = 10"),
  make_df(mgam$results$peak.t, mgam20$results$peak.t, "t", "k = 20"),
  make_df(mgam$results$peak.t, mgam50$results$peak.t, "t", "k = 50"),
  make_df(mgam$results$peak.t, mgam_hd$results$peak.t, "t", "Hand drawn"),

  make_df(mgam$results$peak.r, mgam10$results$peak.r, "r", "k = 10"),
  make_df(mgam$results$peak.r, mgam20$results$peak.r, "r", "k = 20"),
  make_df(mgam$results$peak.r, mgam50$results$peak.r, "r", "k = 50"),
  make_df(mgam$results$peak.r, mgam_hd$results$peak.r, "r", "Hand drawn")
)

df$stat <- factor(df$stat, levels = c("t", "r"))
df$method <- factor(df$method, levels = c("k = 10", "k = 20", "k = 50", "Hand drawn"))

cors <- do.call(
  rbind,
  lapply(split(df, list(df$stat, df$method), drop = TRUE), function(z) {
    data.frame(
      stat = unique(z$stat),
      method = unique(z$method),
      r = cor(z$reference, z$comparison, use = "complete.obs"),
      x = min(z$reference, na.rm = TRUE),
      y = max(z$comparison, na.rm = TRUE)
    )
  })
)

p <- ggplot(df, aes(reference, comparison)) +
  geom_point(size = 1.4, alpha = 0.65) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = cors,
    aes(
      x = x,
      y = y,
      label = paste0("r = ", round(r, 3))
    ),
    hjust = 0,
    vjust = 1,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  facet_grid(
    stat ~ method,
    labeller = labeller(
      stat = c(t = "Peak t", r = "Peak r")
    )
  ) +
  labs(
    x = "Peak statistics with k = 14",
    y = "Peak statistics with other k"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11),
    aspect.ratio = 1
  )



ggsave(p, filename="../plots/mucosa_k_robustness.png",
       width= 6.5,
       height= 6.5)