
library(MorphoGAM)
library(RColorBrewer)

setwd(here::here("analysis/moffitt_mucosa/code"))

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























library(ggrepel)
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

ggsave(p.peak, filename="../plots/mouse_mucosa_peak.png")

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

fpc1 <- plotFPCloading(mgam_object=mgam,
                       curve.fit=fit,
                       L=1,num_genes=5)

fpc2 <- plotFPCloading(mgam_object=mgam,
                       curve.fit=fit,
                       L=2,num_genes=4)




df.full <- readRDS("../data/gsea_results.RDS")


set.seed(1)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggtext)

# Step 1: Define always-include gene sets
always_include <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_APICAL_SURFACE",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_PANCREAS_BETA_CELLS",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
)

# Step 2: Top 2 gene sets by statistic (x) per method
top2_per_method <- df.full %>%
  group_by(Method) %>%
  slice_max(order_by = x, n = 2) %>%
  ungroup() %>%
  pull(y)

# Step 3: Filter to keep only relevant gene sets
keep_genes <- union(always_include, top2_per_method)

df.filtered <- df.full %>%
  filter(y %in% keep_genes)

# Step 4: Set ordering and HTML label coloring
highlighted <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")
set.seed(123)
other_genes <- setdiff(unique(df.filtered$y), highlighted)
ordered_genes <- c(highlighted, sample(other_genes))

df.filtered <- df.filtered %>%
  mutate(
    y = factor(y, levels = ordered_genes),
    # Remove 'HALLMARK_' from display label
    clean_label = gsub("^HALLMARK_", "", as.character(y)),
    label = ifelse(as.character(y) %in% highlighted,
                   paste0("<span style='color:red'>", clean_label, "</span>"),
                   clean_label)
  )

p.sub <- ggplot(data = df.filtered, aes(x = x, y = fct_rev(label))) +
  geom_segment(aes(x = 0, xend = x, y = fct_rev(label), yend = fct_rev(label)),
               color = "gray60", linewidth = 0.5) +
  geom_point(aes(color = y %in% highlighted), size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
  facet_wrap(~Method, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_markdown(),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "GSEA Test Statistic", y = NULL)








library(ggpubr)
p <- ggarrange(p.peak, p.range, nrow=2, labels=c("a","b"))

p2 <- ggarrange(p, p.sub, ggarrange(fpc1,fpc2, nrow=1), nrow=3,
                heights=c(2,1,1), labels=c("", "c", "d"))


ggsave(p2, filename="../plots/mucosa_updated_gsea_and_pcs.png",
       width= 4.78*1.5,
       height= 8.18*1.5,
       units="in")


