setwd(here::here("analysis/moffitt_mucosa/code"))

library(tidyverse)

library(piano)

#gene_sets <- loadGSC(file="../data/m5.go.bp.v2025.1.Mm.symbols.gmt",
#                     type="gmt")

gene_sets <- loadGSC(file="../data/mh.all.v2025.1.Mm.symbols.gmt",
                     type="gmt")


# View the top enriched Hallmark pathways
#print(head(summary_table_sorted, 10))

nnsvg <- readRDS("../data/nnSVG_results.RDS")

#nnsvg <- nnsvg[1:100,]
ranked.list <- nnsvg$LR_stat
names(ranked.list) <- rownames(nnsvg)
#names(ranked.list) <- toupper(rownames(nnsvg))


gsa_result <- runGSA(
  geneLevelStats = ranked.list,
  geneSetStat = "mean",       # or "mean", "median", "fisher", etc.
  gsc = gene_sets,
  nPerm = 1000                # Increase for stability
)

print(GSAsummaryTable(gsa_result))

summary_table <- GSAsummaryTable(gsa_result)

summary_table_sorted <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

# View the top enriched Hallmark pathways
print(head(summary_table_sorted, 10))

summary_table <- GSAsummaryTable(gsa_result)

gsea.nnsvg <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

# View the top enriched Hallmark pathways
#print(head(summary_table_sorted, 10))

df.nnsvg <- data.frame(y=gsea.nnsvg$Name,
                      x=gsea.nnsvg$`Stat (non-dir.)`,
                      Method="nnSVG")



spark <- readRDS("../data/spark_results.RDS")



#spark <- spark[1:100,]
ranked.list <- -log10(spark$combinedPval)
#names(ranked.list) <- toupper(rownames(spark))
names(ranked.list) <- rownames(spark)

gsa_result <- runGSA(
  geneLevelStats = ranked.list,
  geneSetStat = "mean",       # or "mean", "median", "fisher", etc.
  gsc = gene_sets,
  nPerm = 1000                # Increase for stability
)

print(GSAsummaryTable(gsa_result))

summary_table <- GSAsummaryTable(gsa_result)

summary_table_sorted <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

# View the top enriched Hallmark pathways
print(head(summary_table_sorted, 10))

summary_table <- GSAsummaryTable(gsa_result)

gsea.spark <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

# View the top enriched Hallmark pathways
#print(head(summary_table_sorted, 10))


df.spark <- data.frame(y=gsea.spark$Name,
                       x=gsea.spark$`Stat (non-dir.)`,
                       Method="Spark-X")



load("../data/mucosa_mgam.RData")

mgam.df <- mgam$results |> arrange(desc(peak.t))

#mgam.df <- mgam.df[1:100,]

ranked.list <- mgam.df$peak.t
#names(ranked.list) <- toupper(rownames(mgam.df))
names(ranked.list) <- rownames(mgam.df)


gsa_result <- runGSA(
  geneLevelStats = ranked.list,
  geneSetStat = "mean",       # or "mean", "median", "fisher", etc.
  gsc = gene_sets,
  nPerm = 1000                # Increase for stability
)

summary_table <- GSAsummaryTable(gsa_result)

gsea.morphogam <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]


df.mgam <- data.frame(y=gsea.morphogam$Name,
                      x=gsea.morphogam$`Stat (non-dir.)`,
                      Method="MorphoGAM")


df.full <- rbind(df.spark, df.mgam, df.nnsvg)

saveRDS(df.full, "../data/gsea_results.RDS")

p <- ggplot(data=df.full, aes(x=x,y=y)) +
  geom_point() + facet_wrap(~Method, scales="free_x")


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
    label = ifelse(as.character(y) %in% highlighted,
                   paste0("<span style='color:red'>", as.character(y), "</span>"),
                   as.character(y))
  )

# Step 5: Plot with lollipops
p.sub <- ggplot(data = df.filtered, aes(x = x, y = fct_rev(label))) +
  # Lollipop lines
  geom_segment(aes(x = 0, xend = x, y = fct_rev(label), yend = fct_rev(label)),
               color = "gray60", linewidth = 0.5) +
  # Points
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
  labs(x = "GSEA TODO", y = NULL)




library(ggplot2)
library(dplyr)
library(forcats)
library(ggtext)

# Step 1: Use the full dataset (no filtering)
df.all <- df.full

# Step 2: Highlight the same interferon gene sets
highlighted <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")

# Step 3: Set order: interferons on top, rest randomized
set.seed(123)
all_genes <- unique(df.all$y)
other_genes <- setdiff(all_genes, highlighted)
ordered_genes <- c(highlighted, sample(other_genes))

df.all <- df.all %>%
  mutate(
    y = factor(y, levels = ordered_genes),
    label = ifelse(as.character(y) %in% highlighted,
                   paste0("<span style='color:red'>", as.character(y), "</span>"),
                   as.character(y))
  )

# Step 4: Plot with lollipops for all gene sets
p.full <- ggplot(data = df.all, aes(x = x, y = fct_rev(label))) +
  # Lollipop lines
  geom_segment(aes(x = 0, xend = x, y = fct_rev(label), yend = fct_rev(label)),
               color = "gray60", linewidth = 0.5) +
  # Points
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
  labs(x = "GSEA test statistic", y = NULL)


ggsave(p.full,
       filename="../plots/gsea_full.png",
       width=1.5*6.73, height=1.5*3.87)

