

library(STexampleData)

spe <- SlideSeqV2_mouseHPC()


ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]

out <- CurveSearcher(xy,knn=5)

out$plot

ggsave(out$plot, file="../plots/curve.png")

save(out, file="../data/curve_object.rda")
