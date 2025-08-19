setwd(here::here("analysis/gaston/code"))

##Run moffitt mucosa to get Y.sub


isodepth <- read.csv(file="../data/gaston_isodepth.csv")

isodepth <- as.vector(isodepth[,1])

Y.sub <- Y.sub[,-1]





plot(isodepth, Y.sub["Apob",])

plot(isodepth, Y.sub["Ddx58",])



library(tidyverse)
library(ggpubr)

locus <- as.matrix(meta.sub[,c("x","y")])

#fit <- CurveFinder(locus,knn=10,loop=TRUE)

gene.1 <- Y.sub["Ddx58",]
gene.2 <- Y.sub["Apob",]

p.gene1 <- data.frame(x=isodepth, y=gene.1) |>
  ggplot(aes(x=x,y=y)) +
  geom_jitter(width=0,height=0.1,size=0.5) +
  theme_bw() +
  xlab("Isodepth") + ylab("Count") + ggtitle("Ddx58")

p.gene2 <- data.frame(x=isodepth, y=gene.2) |>
  ggplot(aes(x=x,y=y)) +
  geom_jitter(width=0,height=0.1,size=0.5) +
  theme_bw() +
  xlim(0.3,1) +
  xlab("Isodepth") + ylab("Count") + ggtitle("Apob")

p <- ggarrange(p.gene1, p.gene2)

ggsave(p, filename="../plots/isodepth_gene_plots.png",
       width=10.6, height=5.6, units="in")

