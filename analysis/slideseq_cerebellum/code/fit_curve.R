
data.dir <- "/Volumes/pnicol_data/MERFISH/SCP948"
res <- readRDS(file.path(data.dir, "seurat/myRCTD_cerebellum_slideseq.rds"))
df <- res@results$results_df
xy <- readRDS(file.path(data.dir, "seurat/puckCropped_cerebellum_slideseq.rds"))
xy <- xy@coords
df$barcode <- rownames(df)
xy$barcode <- rownames(xy)

xy <- left_join(xy, df, by="barcode")
xy <- xy[-which(is.na(xy$first_type)),]
xy <- xy[xy$first_type == "Granule", ]
xy <- xy[,c("x","y")]
xy <- as.matrix(xy)

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
       file="../plots/ca3_cerebellum.png")

save(out, file="../data/curve_object.rda")
