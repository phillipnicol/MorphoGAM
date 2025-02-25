setwd(here::here("analysis/moffitt_mucosa/code"))

Y <- read.csv("../../data/GL2_distal_colon_cell_by_gene_raw.csv")

rownames(Y) <- Y[,1]

Y <- Y[,-1]

library(Matrix)

Y <- as.matrix(Y)

Y <- t(Y) #Transpose to get genes x cells


meta <- read.csv("../../data/GL2_distal_colon_cell_type_and_locations_2023.08.11.csv")

library(tidyverse)
meta <- as.data.frame(meta)
meta.sub <- meta |> filter(slice_full_name == "20220518_WT_dcol_slice_3") |>
  filter(spatial_neighborhood_v1 == "Mucosa") |>
  filter(leiden_combined_v2 == "Enterocyte")

Y.sub <- Y[,meta.sub$X]



gene.1 <- which(rownames(Y.sub) == "Apob")
gene.2 <- which(rownames(Y.sub) == "Ddx58")

#Plot these two
#expr <- t(Y.sub[c(gene.1,gene.2),]) |>
# apply(2,function(x) log2(x+1)) |>
#  as.data.frame() |>
#  mutate(x=meta.sub$x, y=meta.sub$y) |>
#  pivot_longer(cols=-c(x,y))

#plot1 <- expr |> ggplot(aes(x=x,y=y,color=value)) +
#  geom_point(size=0.25,alpha=0.75) +
#  scale_color_gradient(low="grey85", high="darkred")+
#  facet_wrap(~name) +
#  labs(color="log expression") +
#  theme_bw()

Y.s <- sweep(Y.sub, MARGIN = 2, STATS = colSums(Y.sub)/median(colSums(Y.sub)),
             FUN="/")

#Plot these two
expr <- t(Y.s[c(gene.1,gene.2),]) |>
  apply(2,function(x) log2(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))
expr$name <- factor(expr$name, levels=c("Ddx58", "Apob"))

plot1 <- expr |> ggplot(aes(x=x,y=y,color=value,size=value,
                            alpha=value)) +
  geom_point(size=0.25) +
  scale_color_gradient(low="grey85", high="darkred")+
  #scale_size_continuous(range=c(0.1,0.75)) +
  scale_alpha_continuous(range=c(0.33,1)) +
  coord_fixed() +
  facet_wrap(~name) +
  labs(color="log expression") +
  theme_bw() +
  guides(alpha="none")


#ggsave(plot1,"../plots/two_genes_examples.png")




gene.1 <- which(rownames(Y.sub) == "Tnfsf10")
gene.2 <- which(rownames(Y.sub) == "Nlrc5")
#Plot these two
expr <- t(Y.sub[c(gene.1,gene.2),]) |>
  apply(2,function(x) log2(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y)) |>
  mutate(name=ifelse(name == "Tnfsf10",
                     "Tnsfs10
                     (Spark Rank X),
                     (nnSVG Rank X)",
                     "Nlrc5
                     (Spark Rank X),
                     (nnSVG rank X)"))

plot2 <- expr |>
  ggplot(aes(x=x, y=y, color=value)) +
  # Layer for zero values
  geom_point(data = ~ subset(., value == 0),
             size = 1, alpha = 0.5) +
  # Layer for non-zero values
  geom_point(data = ~ subset(., value != 0),
             size = 1.5, alpha = 1) +
  scale_color_gradient(low="grey90", high="darkred") +
  facet_wrap(~name) +
  xlim(c(9000, 9500)) + ylim(c(-750,250)) +
  labs(color="log expression") +
  theme_bw()

ggsave(plot2,filename="../plots/two_other_examples.png")




### SPARKX

#Test expression for these genes
library(SPARK)
locus <- as.matrix(meta.sub[,c("x","y")])
res <- SPARK::sparkx(count_in = Y.sub,
                     locus_in = locus)

p.vals <- res$res_mtest |> as.data.frame() |>
  arrange(adjustedPval)

## Take top 12
top.12 <- rownames(p.vals)[1:12]

#Plot these two
expr <- t(Y.sub[top.12,]) |>
  apply(2,function(x) log2(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))

plot2 <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.25,alpha=0.75) +
  scale_color_gradient(low="grey90", high="darkred")+
  coord_fixed() +
  facet_wrap(~name,nrow=4,ncol=3) +
  labs(color="log expression") +
  theme_bw()
ggsave(plot2,filename="../plots/spark_top12.png",width=8,height=10)

print(res$res_mtest["Ddx58",])
print(res$res_mtest["Apob",])


#logCPM <- log(sweep(Y.sub, MARGIN=2, STATS = colSums(Y.sub), FUN="/")+1)
#nn.svg <- nnSVG(input=logCPM, spatial_coords = locus,
#                verbose=TRUE)

nn.svg <- readRDS("../data/nnSVG_results.RDS")

p.vals.nnsvg <- nn.svg |> as.data.frame() |>
  arrange(rank)

top.12.nnsvg <- rownames(p.vals.nnsvg)[1:12]
expr <- t(Y.sub[top.12.nnsvg,]) |>
  apply(2,function(x) log2(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))
plot3 <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.25,alpha=0.75) +
  scale_color_gradient(low="grey90", high="darkred")+
  facet_wrap(~name,nrow=4,ncol=3) +
  labs(color="log expression") +
  theme_bw()

print(nn.svg["Ddx58",])
print(nn.svg["Apob",])

ggsave(plot3,filename="../plots/nnSVG_top12.png",width=8,height=10)






## Other examples of 1D paths
p.mucosa <- expr |> ggplot(aes(x=x,y=y)) +
  geom_point(size=0.5) +
  theme_bw() +
  ggtitle("MERFISH mouse mucosa enterocytes") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") + ylab("")



## Granule

xy <- readRDS("../../data/granule_slideseq.RDS")

## Prune outlier
xy.dist <- as.matrix(dist(xy))
nnk <- apply(xy.dist, 1, function(x) sort(x)[6])
outlier <- which(nnk > 2*median(nnk))
xy <- xy[-outlier,]


p.granule <- xy |> ggplot(aes(x=x,y=y)) + geom_point(size=0.5) +
  theme_bw() +
  ggtitle("Slide-seq cerebellum granule cells") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  theme(plot.title = element_text(size = 11))



## CA3

spe <- STexampleData::SlideSeqV2_mouseHPC()

ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]
Y <- counts(spe)[,ixs]

xy.dist <- as.matrix(dist(xy))
knn <- 20
prune.outlier <- 3
nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
outlier <- which(nnk > (prune.outlier)*median(nnk))

xy <- xy[-outlier,]

p.ca3 <- xy |> ggplot(aes(x=xcoord,y=ycoord)) + geom_point(size=0.5) +
  theme_bw() +
  ggtitle("Slide-seq hippocampus CA3 cells") +
  xlab("") + ylab("") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  theme(plot.title = element_text(size = 11))

library(ggpubr)

plot1 <- plot1 + ggtitle("MERFISH mouse colon enterocytes") +
  theme(plot.title = element_text(size = 11))

### Also add swiss roll examples here


load("../../merfish_swissroll/data/curve_d9_061923.RData")

p.swiss1 <- data.frame(x=fit$xyt$x, y=fit$xyt$y) |>
  ggplot(aes(x=x,y=y)) + geom_point(size=0.25) +
  theme_bw() +
  ggtitle("MERFISH swissroll 1") +
  xlab("") + ylab("") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  theme(plot.title = element_text(size = 11))



load("../../merfish_swissroll/data/curve_d9_m5_080823.RData")

p.swiss2 <- data.frame(x=fit$xyt$x, y=fit$xyt$y) |>
  ggplot(aes(x=x,y=y)) + geom_point(size=0.25) +
  theme_bw() +
  ggtitle("MERFISH swissroll 2") +
  xlab("") + ylab("") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  theme(plot.title = element_text(size = 11))

load("../../merfish_swissroll/data/curve_d9_m13_080823.RData")

p.swiss3 <- data.frame(x=fit$xyt$x, y=fit$xyt$y) |>
  ggplot(aes(x=x,y=y)) + geom_point(size=0.25) +
  theme_bw() +
  ggtitle("MERFISH swissroll 3") +
  xlab("") + ylab("") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  theme(plot.title = element_text(size = 11))

p <- ggarrange(plot1, ggarrange(p.granule, p.ca3,nrow=1, labels=c("b","c")),
               ggarrange(p.swiss1, p.swiss2, p.swiss3,nrow=1),
               nrow=3, labels=c("a","", "d"),
               heights=c(1.5,1,1))


ggsave(p, filename="../plots/example_of_1d_new.png",
       width=7.56, height=9.75, units="in")

#ggsave(p, filename="../plots/example_of_1d.png",
#       width=11.88, height=9.19, units="in")



## Curve Demo

library(MorphoGAM)
library(tidyverse)
library(ggpubr)

locus <- as.matrix(meta.sub[,c("x","y")])

fit <- CurveFinder(locus,knn=10,loop=TRUE)

gene.1 <- Y.sub["Ddx58",]
gene.2 <- Y.sub["Apob",]

p.gene1 <- data.frame(x=fit$xyt$t, y=gene.1) |>
  ggplot(aes(x=x,y=y)) +
  geom_jitter(width=0,height=0.1,size=0.5) +
  theme_bw() +
  xlab("t") + ylab("Count") + ggtitle("Ddx58")

p.gene1r <- data.frame(x=fit$xyt$r, y=gene.1) |>
  ggplot(aes(x=x,y=y)) +
  geom_jitter(width=0,height=0.1,size=0.5) +
  theme_bw() +
  xlim(0.3,1) +
  xlab("r") + ylab("Count") + ggtitle("Ddx58")

p.gene2 <- data.frame(x=fit$xyt$r, y=gene.2) |>
  ggplot(aes(x=x,y=y)) +
  geom_jitter(width=0,height=0.1,size=0.5) +
  theme_bw() +
  xlim(0.3,1) +
  xlab("r") + ylab("Count") + ggtitle("Apob")


p.gene2t <- data.frame(x=fit$xyt$t, y=gene.2) |>
  ggplot(aes(x=x,y=y)) +
  geom_jitter(width=0,height=0.1,size=0.5) +
  theme_bw() +
  xlab("t") + ylab("Count") + ggtitle("Apob")

p.demo <- ggarrange(ggarrange(fit$curve.plot + ggtitle("Fitted curve") + guides(color="none"),
                              fit$coordinate.plot + guides(color="none") + ggtitle("First coordinate (t)"),
                              fit$residuals.plot + guides(color="none") + ggtitle("Second coordinate (r)"), nrow=1,ncol=3,
                              labels=c("a","b","c")),
                    ggarrange(p.gene1, p.gene2,
                              p.gene1r, p.gene2t,
                              nrow=2,ncol=2,
                              labels=c("d","e", "f","g")),nrow=2,
                    heights=c(1,2))

ggsave(p.demo, filename="../plots/curve_demo.png",
       width=8.89, height=9.54, units="in")


