library(tidyverse)
meta <- read.csv("../../data/GL2_distal_colon_cell_type_and_locations_2023.08.11.csv")
meta <- as.data.frame(meta)
meta.sub <- meta |> filter(slice_full_name == "20220518_WT_dcol_slice_3") |>
  filter(spatial_neighborhood_v1 == "Mucosa") |>
  filter(leiden_combined_v2 == "Enterocyte")

xy <- as.matrix(meta.sub[,c("x","y")])

library(MorphoGAM)

#fit <- CurveFinderInteractive(xy,loop=TRUE)
#save(fit, file="../data/loop_ground_truth.Rda")
load("../data/loop_ground_truth.Rda")



p <- xy |> ggplot(aes(x=x,y=y)) + geom_point(size=0.5) +
  theme_bw() +
  ggtitle("Mouse colon")
p1 <- fit$curve.plot + guides(color="none") + ggtitle("Path (Hand drawn)")
p2 <- fit$coordinate.plot +guides(color="none")


out10 <- CurveFinder(xy,knn=10, loop=TRUE)
p.cf.10 <- out10$curve.plot + guides(color="none") + ggtitle("MorphoGAM (knn=10)")

out15 <- CurveFinder(xy,knn=10, loop=TRUE)
p.cf.15 <- out15$curve.plot + guides(color="none")+ ggtitle("MorphoGAM (knn=15)")

out20 <- CurveFinder(xy,knn=10, loop=TRUE)
p.cf.20 <- out20$curve.plot + guides(color="none")+ ggtitle("MorphoGAM (knn=20)")


### Principal curves

library(princurve)

fit.ps <- principal_curve(as.matrix(xy), smoother="periodic_lowess")

fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)
fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p.ps.1 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("princurve (Default)") +
  guides(color="none")

fit.ps <- principal_curve(as.matrix(xy),smoother="periodic_lowess",
                          f=1/3)

fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)
fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p.ps.2 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("princurve (f=1/3)") +
  guides(color="none")

fit.ps <- principal_curve(as.matrix(xy),smoother="periodic_lowess",
                          f=0.01)

fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)
fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p.ps.3 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("princurve (f=0.01)") +
  guides(color="none")

library(ggpubr)

p.big <- ggarrange(ggarrange(p,p1,p2,nrow=1,ncol=3),
                   ggarrange(p.cf.10, p.cf.15, p.cf.20,nrow=1,ncol=3),
                   ggarrange(p.ps.1,p.ps.2,p.ps.3, nrow=1,ncol=3),
                   nrow=3)

ggsave(p.big, filename="../plots/loop.png",width=8, height=8, units="in")
