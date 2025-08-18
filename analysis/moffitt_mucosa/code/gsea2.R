setwd(here::here("analysis/moffitt_mucosa/code"))

library(tidyverse)

library(piano)

gene_sets <- loadGSC(file="../data/m5.go.bp.v2025.1.Mm.symbols.gmt",
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
                       x=-log10(gsea.nnsvg$`p (non-dir.)`),
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
                       x=-log10(gsea.spark$`p (non-dir.)`),
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
                      x=-log10(gsea.morphogam$`p (non-dir.)`),
                      Method="MorphoGAM")


df.full <- rbind(df.spark, df.mgam, df.nnsvg)

saveRDS(df.full, "../data/gsea2_results.RDS")




library(dplyr)
library(ggplot2)
library(stringr)

# Mark interferon-related GO terms (tweak the pattern if you want it stricter/looser)
df_plot <- df.full %>%
  mutate(
    Interferon = str_detect(
      y,
      regex("interferon|\\bIFN\\b|type\\s*[I1]\\s*interferon|type\\s*II\\s*interferon",
            ignore_case = TRUE)
    )
  )

# Optional: order methods the way you prefer
# df_plot <- df_plot %>%
#   mutate(Method = factor(Method, levels = c("Spark-X", "MorphoGAM", "nnSVG")))

ggplot(df_plot, aes(x = Method, y = x)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey40") +
  geom_jitter(
    data = dplyr::filter(df_plot, Interferon),
    width = 0.18, height = 0,
    color = "red", size = 2.5
  ) +
  labs(
    x = NULL,
    y = "-log10 p-value",
    title = "GO Biological Process Enrichment by Method",
    subtitle = "Red points = Interferon-related gene sets"
  ) +
  theme_classic(base_size = 12)



library(dplyr)
library(ggplot2)
library(stringr)

df_plot <- df.full %>%
  mutate(
    Interferon = if_else(
      str_detect(
        y,
        regex("interferon|\\bIFN\\b|type\\s*[I1]\\s*interferon|type\\s*II\\s*interferon",
              ignore_case = TRUE)
      ),
      "Interferon gene sets",
      "Other gene sets"
    )
  )

ggplot(df_plot, aes(x = Method, y = x, fill = Interferon)) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "grey30"
  ) +
  geom_jitter(
    color = "black",
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
    alpha = 0.6, size = 1.8
  ) +
  labs(
    x = NULL,
    y = "-log10(p-val)",
    title = "GO Biological Process Enrichment by Method",
    subtitle = "Dodged boxplots for Interferon vs. Other gene sets"
  ) +
  scale_fill_manual(
    values = c("Interferon gene sets" = "#1f78b4",   # blue
               "Other gene sets"      = "grey80")
  ) +
  theme_bw()
