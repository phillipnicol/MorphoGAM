setwd(here::here("analysis/moffitt_mucosa/code"))

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

summary_table_sorted <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

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

summary_table_sorted <- summary_table[order(summary_table$`Stat (non-dir.)`,decreasing = TRUE), ]

# View the top enriched Hallmark pathways
print(head(summary_table_sorted, 10))

