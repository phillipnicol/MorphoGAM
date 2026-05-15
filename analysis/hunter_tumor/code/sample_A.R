setwd(here::here("analysis", "hunter_tumor", "code"))

library(ggplot2)

xy_raw <- read.csv("../data/coords_sampleA.csv", header = FALSE)

full_w <- 31870
full_h <- 32588

small_w <- 785
small_h <- 768

image_path <- "../data/GSM4838131_Visium_A_image.tif"

ixs <- which(xy_raw$V2 == 1)

x_full <- xy_raw[ixs, 5]
y_full <- xy_raw[ixs, 6]

# aspect-preserving shrink matched to height
s <- small_h / full_h

expected_w <- full_w * s
pad_x <- (small_w - expected_w) / 2

# manual correction: dots looked too far right
x_shift <- 0
y_shift <- 10

xy_small <- cbind(
  x = x_full * s + pad_x + x_shift,
  y = y_full * s + y_shift
)

# Do NOT flip y here, because CurveFinderInteractive2 does it
fit <- CurveFinderInteractive2(
  xy = xy_small,
  loop = TRUE,
  image_path = image_path,
  flip_y = FALSE
)





p1 <- fit$curve.plot + ggtitle("Fitted curve along tumor boundary")

t.smooth <- mgcv::gam(fit$xyt$t ~ s(fit$xyt$x, fit$xyt$y, bs="tp", k=100))

fit$xyt$t <- predict(t.smooth)

p2 <- data.frame(x=fit$xyt$x,y=fit$xyt$y,color=fit$xyt$t) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),
                        colors=c("blue","grey90", "red"))+
  theme_bw() +
  ggtitle("Coordinate along curve") + labs(color="t")

r <- ifelse(fit$xyt$r < 0.5, fit$xyt$r - 0.5, -1)

p2 <- data.frame(x=fit$xyt$x,y=fit$xyt$y,color=r) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),
                        colors=c("blue","grey90", "red"))+
  theme_bw() +
  ggtitle("Distance from curve") + labs(color="Distance")





library(ggplot2)

df <- data.frame(x=fit$xyt$x, y=-fit$xyt$y, r=fit$xyt$r)

p <- ggplot() +

  # Tumor (or graded) cells: r < 0.5
  geom_point(
    data = df[df$r > 0.5, ],
    aes(x = x, y = y, color = 0.5 - r),
    size = 1
  ) +

  # Non-tumor cells: r > 0.5
  geom_point(
    data = df[df$r < 0.5, ],
    aes(x = x, y = y),
    color = "grey90",
    size = 1
  ) +

  scale_color_gradient(
    low = "darkblue",
    high = "lightblue",
    name = "r"
  ) +

  guides(
    color = guide_colorbar(order = 1),
    fill  = "none"
  ) +

  labs(
    color = "r",
    caption = "Grey points: non-tumor cell"
  ) +

  theme_bw()



p1 <- fit$curve.plot + ggtitle("Fitted curve along tumor boundary")
p <- p + ggtitle("Second coordinate (distance from tumor edge)")

kxs <- which(fit$xyt$r > 0.5)


library(anndata)

library(rhdf5)

library(rhdf5)
library(Matrix)

file <- "../data/GSM4838132_Visium_Sample_B_filtered_feature_bc_matrix.h5"

data    <- h5read(file, "matrix/data")
indices <- h5read(file, "matrix/indices")
indptr  <- h5read(file, "matrix/indptr")
shape   <- h5read(file, "matrix/shape")

# Build sparse matrix (note: +1 because R is 1-based)
mat <- sparseMatrix(
  i = indices + 1,
  p = indptr,
  x = data,
  dims = shape
)

genes <- h5read(file, "matrix/features/name")
cells <- h5read(file, "matrix/barcodes")

rownames(mat) <- genes
colnames(mat) <- cells


cell_meta <- data.frame(
  barcode = h5read(file, "matrix/barcodes"),
  stringsAsFactors = FALSE
)

head(cell_meta)


mat_r <- as(mat, "RsparseMatrix")

jxs <- match(xy_raw[ixs,]$V1, colnames(mat_r))

mat_r <- mat_r[,jxs]

gene <- mat_r[which(rownames(mat_r) == "RPL41"),]

fit$xyt <- fit$xyt[kxs,]
mat_r <- mat_r[,kxs]


data.frame(x=xy_small[,1], y=xy_small[,2], color=gene) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void()

# Find top 2000 spatially variable genes
#Change this later 
sample.var <- apply(mat_r, 1, var)
top_genes <- names(sort(sample.var, decreasing = TRUE))[1:2000]
mat_r <- mat_r[top_genes,]
rownames(mat_r) <- make.unique(rownames(mat_r))

mgam <- MorphoGAM(Y=mat_r,
                  curve.fit = fit,
                  design = y ~ s(t) + s(r))


plot(fit$xyt$r, mat_r[4,])

peak_r_location <- fit$xyt$r[apply(mgam$fxs.r, 1, which.max)]

rownames(mat_r)[order(abs(peak_r_location - 0.5), decreasing = FALSE)[1:10]]

plot(fit$xyt$r, mat_r["lygl1",])

plot(fit$xyt$t, mat_r["vmhc",])



data.frame(x=fit$xyt$x, y=fit$xyt$y, color=mat_r["pvalb4",]) |>
  ggplot(aes(x=x,y=y,color=color)) + geom_point() +
  scale_color_gradient(low="lightgrey", high="red") + theme_void()


plot(fit$xyt$t, mat_r["pvalb4",])
