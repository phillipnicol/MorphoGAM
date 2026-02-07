setwd(here::here("analysis", "hunter_tumor", "code"))

### SAMPLE B

xy <- read.csv("../data/coords_sampleC.csv", header=F)

#31439x32445
full_w <- 31439
full_h <- 32445

#785x768
small_w <- 792
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
  loop = TRUE,
  image_path = "/Users/phillipnicol/Desktop/GSM4838132_Visium_B_image.tif",
  flip_y = TRUE
)


p1 <- fit$curve.plot + ggtitle("Fitted curve along tumor boundary")

t.smooth <- mgcv::gam(fit$xyt$t ~ s(fit$xyt$x, fit$xyt$y, bs="tp", k=100))

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
    data = df[df$r < 0.5, ],
    aes(x = x, y = y, color = 0.5 - r),
    size = 1
  ) +

  # Non-tumor cells: r > 0.5
  geom_point(
    data = df[df$r > 0.5, ],
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

library(ggpubr)

p.big <- ggarrange(p1, p)

ggsave(p.big, filename="../plots/sampleB.png",
       width=9.25, height=4.41, units="in")


p1

ggsave(p1, filename="../plots/sample_b_curve.png")
ggsave(p, filename="../plots/sample_b_distance.png")

