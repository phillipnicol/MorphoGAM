setwd(here::here("analysis/moffitt_mucosa/code"))

library(tidyverse)

library(piano)

load("../data/mucosa_mgam.RData")

mgam.df <- mgam$results |> arrange(desc(peak.t))

mgam.df <- mgam.df[1:100,]

ranked.list <- mgam.df$peak.t
names(ranked.list) <- toupper(rownames(mgam.df))

gene_sets <- loadGSC(file="../data/h.all.v2025.1.Hs.symbols.gmt",
                     type="gmt")

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


# View the top enriched Hallmark pathways
print(head(summary_table_sorted, 10))

nnsvg <- readRDS("../data/nnSVG_results.RDS")

nnsvg <- nnsvg[1:100,]
ranked.list <- nnsvg$LR_stat
names(ranked.list) <- toupper(rownames(nnsvg))


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



spark <- spark[1:100,]
ranked.list <- spark$combinedPval
names(ranked.list) <- toupper(rownames(spark))


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


df.full <- rbind(df.spark, df.mgam, df.nnsvg)

p <- ggplot(data=df.full, aes(x=x,y=y)) +
  geom_point() + facet_wrap(~Method, scales="free_x")


set.seed(1)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggtext)

# Step 1: Define highlighted gene sets
highlighted <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")

# Step 2: Randomize other gene sets
set.seed(123)
all_genes <- unique(df.full$y)
other_genes <- setdiff(all_genes, highlighted)
shuffled_others <- sample(other_genes)

# Step 3: Final ordering: highlighted at the top
ordered_genes <- c(highlighted, shuffled_others)

# Step 4: Apply reordering and HTML labels
df.full <- df.full %>%
  mutate(
    y = factor(y, levels = ordered_genes),  # NO rev() here to keep top at top
    label = ifelse(as.character(y) %in% highlighted,
                   paste0("<span style='color:red'>", as.character(y), "</span>"),
                   as.character(y))
  )

# Step 5: Plot
ggplot(data = df.full, aes(x = x, y = fct_rev(label))) +  # flip y-axis to show top-to-bottom
  geom_point(aes(color = y %in% highlighted)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
  facet_wrap(~Method, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_markdown(),  # enables red text rendering
    strip.text = element_text(face = "bold")
  ) +
  labs(x = NULL, y = NULL)
