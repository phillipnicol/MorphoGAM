

library(STexampleData)


spe <- Visium_humanDLPFC()

for(i in 1:6) {
  print(i)
  layer <- paste0("Layer", i)
  ixs <- which(spe$ground_truth == layer) #subset to Layer

  xy <- spatialCoords(spe)[ixs,]

  out <- CurveSearcher(xy,knn=5)

  out$plot

  ggsave(out$plot, file=paste0("../plots/curve_", layer, ".png"))

  save(out, file=paste0("../data/curve_object_", layer, ".rda"))
}
