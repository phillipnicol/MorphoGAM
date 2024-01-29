
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

out <- CurveSearcher(xy,knn=5)

out$plot

ggsave(out$plot, file="../plots/curve.png")

save(out, file="../data/curve_object.rda")
