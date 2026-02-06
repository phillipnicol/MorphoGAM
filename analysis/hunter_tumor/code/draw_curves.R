setwd(here::here("analysis", "hunter_tumor", "code"))

### SAMPLE A

xy <- read.csv("../data/coords_sampleA.csv", header=F)

#31439x32445
full_w <- 31439
full_h <- 32445

#785x768
small_w <- 785
small_h <- 768

ixs <- which(xy$V2 == 1)

xy <- xy[ixs,c(5,6)]

#info_small <- image_info(image_read("/Users/phillipnicol/Desktop/GSM4838131_Visium_A_image.tif"))
#small_w <- info_small$width
#small_h <- info_small$height

sx <- small_w / full_w
sy <- small_h / full_h


x_full <- xy[,1]
y_full <- xy[,2]

xy_small <- cbind(
  x_full * sx,
  y_full * sy
)

xy_small[,2] <- small_h - xy_small[,2]


fit <- CurveFinderInteractive2(
  xy = xy_small,
  loop = FALSE,
  image_path = "/Users/phillipnicol/Desktop/GSM4838131_Visium_A_image.tif",
  flip_y = TRUE
)
