

library(STexampleData)


spe <- Visium_humanDLPFC()

for(i in 1:6) {
  print(i)
  layer <- paste0("Layer", i)
  ixs <- which(spe$ground_truth == layer) #subset to Layer

  xy <- spatialCoords(spe)[ixs,]

  out <- CurveSearcher(xy,knn=5,tau=100)

  p <- out$plot+ggtitle(paste("DLPFC layer ", i)) + guides(color="none")
  p <- p + theme_void()
  p <- p + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  ggsave(p,
         file=paste0("../plots/curve_", layer, ".png"))

  save(out, file=paste0("../data/curve_object_", layer, ".rda"))
}
