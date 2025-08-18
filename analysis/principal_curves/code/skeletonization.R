setwd(here::here("analysis/principal_curves/code"))

# Skeletonize x–y points using {magick}
# pts: data.frame with columns x, y (arbitrary ranges)

library(magick)
library(scales)

#This function was written with chatgpt5
skeletonize_points_magick <- function(pts,
                                      nx = 512, ny = 512,
                                      point_cex = 0.8,
                                      close_disk = 3,
                                      dilate_disk = 2,
                                      thin_iters = 80) {
  stopifnot(all(c("x","y") %in% names(pts)))

  # 1) Normalize to pixel coords
  xs <- rescale(pts$x, to = c(1, nx))
  ys <- rescale(pts$y, to = c(1, ny))

  # 2) Open a magick graphics device and KEEP the returned image object
  img <- image_graph(width = nx, height = ny, bg = "black")
  par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(1, nx), ylim = c(1, ny),
       xlab = "", ylab = "", axes = FALSE)
  points(xs, ys, pch = 16, cex = point_cex, col = "white")
  dev.off()  # finalize (does not return the image)

  # 3) Binarize & clean
  I <- image_convert(img, colorspace = "gray")
  I <- image_threshold(I, type = "white", threshold = "50%")
  I <- image_morphology(I, method = "Close",  kernel = sprintf("Disk:%d", close_disk))
  I <- image_morphology(I, method = "Dilate", kernel = sprintf("Disk:%d", dilate_disk))
  I <- image_threshold(I, type = "white", threshold = "50%")

  # 4) Thin to a 1-px skeleton
  skel <- image_morphology(I, method = "Thinning", kernel = "Skeleton",
                           iterations = thin_iters)
  # Alternative:
  # skel <- image_morphology(I, "Thinning", "ThinSE:8", iterations = thin_iters)

  list(binary = I, skeleton = skel)
}

# --- Example ---
 #pts <- data.frame(x = runif(200), y = runif(200))

## Granule

xy <- readRDS("../../data/granule_slideseq.RDS")

## Prune outlier
xy.dist <- as.matrix(dist(xy))
nnk <- apply(xy.dist, 1, function(x) sort(x)[6])
outlier <- which(nnk > 2*median(nnk))
xy <- xy[-outlier,]
pts <- data.frame(x=xy[,1], y=xy[,2])
 out <- skeletonize_points_magick(pts, nx=512, ny=512)
 print(out$binary); print(out$skeleton)

 image_write(out$skeleton, path = "../plots/granule_skeleton.png", format = "png")

