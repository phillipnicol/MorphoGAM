
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



for(i in 2:6) {
  print(i)
  layer <- paste0("Layer", i)
  ixs <- which(spe$ground_truth == layer) #subset to Layer

  xy <- spatialCoords(spe)[ixs,]

  fit <- CurveFinder(xy, knn=10)

  Y <- Y.big[,ixs]

  mgam <- MorphoGAM(Y = Y,
                    curve.fit = fit,
                    design = y ~ s(t, bs="cr", k=10) + s(r, bs="cr", k=10))

  #my.t <- mgam$results |> arrange(desc(peak.t))

  #plotGAMestimates(Y, genes=rownames(my.t)[1:5], curve_fit=fit, mgam_object = mgam,nrow=1)

  save(fit, file=paste0("../data/curve_fit_layer_", i, ".Rdata"))
  save(mgam, file=paste0("../data/mgam_layer_", i, ".Rdata"))
}
