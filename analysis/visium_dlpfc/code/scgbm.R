

dim(Y.input)

library(scGBM)


scgbm <- gbm.sc(Y.input, M=20)

library(ggpubr)

p.scgbm1 <- data.frame(x=fit$xyt$r, y=scgbm$scores[,1]) |>
  ggplot(aes(x=x,y=y)) + geom_point(size=0.5) +
  theme_bw() + xlab("r") + ylab("PC 1")


p.scgbm2 <- data.frame(x=fit$xyt$r, y=scgbm$scores[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(size=0.5) +
  theme_bw() + xlab("r") + ylab("PC 2")

p.gbmall <- ggarrange(p.scgbm1, p.scgbm2, nrow=1)

#Add title for p.gbmall
p.gbmall <- annotate_figure(p.gbmall, top=text_grob("PCA scores vs r", size=14))

ggsave(p.gbmall, filename="../plots/dlpfc_scgbm.png",
       width=7.23, height=1.89, units="in")



