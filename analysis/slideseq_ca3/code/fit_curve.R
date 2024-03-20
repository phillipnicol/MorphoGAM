

library(STexampleData)

spe <- SlideSeqV2_mouseHPC()


ixs <- which(spe$celltype == "Astrocyte") #subset to CA3

xy <- spatialCoords(spe)[ixs,]

out <- CurveSearcher(xy,knn=5,tau=100)


p <- out$plot+ggtitle("CA3 cells, Mouse HPC") + guides(color="none")
p <- p + theme_void()
p <- p + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
ggsave(p,
       file="../plots/curve_ca3.png")

save(out, file="../data/curve_object.rda")

