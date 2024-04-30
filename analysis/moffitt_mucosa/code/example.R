

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

#Test expression for these genes
library(SPARK)
locus <- as.matrix(meta.sub[,c("x","y")])
res <- SPARK::sparkx(count_in = Y.sub,
                     locus_in = locus)

colors <- c("red", "orange", "yellow", "green", "blue", "purple", "violet")
values <- seq(0,1,length.out=7)
# Plot with custom color scale
ggplot(df, aes(x, y, color = fit$t)) +
  geom_point() +
  scale_color_gradientn(colors = colors,
                        values = values,
                        guide = guide_colorbar(title = "Value",
                                               barheight = 10,
                                               barwidth = 0.5))


df <- data.frame(x=locus[,1],y=locus[,2], color=pst)
ggplot(df, aes(x, y, color = Lineage1)) +
  geom_point() +
  scale_color_gradientn(colors = colors,
                        values = values,
                        guide = guide_colorbar(title = "Value",
                                               barheight = 10,
                                               barwidth = 0.5))

df <- data.frame(x=locus[,1],y=locus[,2], color=fit$t)
ggplot(df, aes(x, y, color = color)) +
  geom_point() +
  scale_color_gradientn(colors = colors,
                        values = values,
                        guide = guide_colorbar(title = "Value",
                                               barheight = 10,
                                               barwidth = 0.5))

p.vals <- res$res_mtest |> as.data.frame() |>
  arrange(adjustedPval)

fit <- CurveSearcherLoop(locus,knn=10)
cso <- fit
my.svg <- detectSVGLoop(Y.sub, cso)

gene.1 <- which(rownames(Y.sub) == "Dhx58")
gene.2 <- which(rownames(Y.sub) == "Erbb3")

#Plot these two
expr <- t(Y.sub[c(gene.1,gene.2),]) |>
  apply(2,function(x) log2(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))

plot1 <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.25,alpha=0.75) +
  scale_color_gradient(low="grey90", high="darkblue")+
  facet_wrap(~name) +
  labs(color="log expression") +
  theme_bw()

##Plot their functions now
expr <- t(Y.sub[c(gene.1,gene.2),-cso$outlier]) |>
  apply(2, function(x) log(x+1)) |>
  as.data.frame() |>
  mutate(x=cso$t) |>
  pivot_longer(cols=-x)


plot3 <- expr |> ggplot(aes(x=x,y=value)) +
  geom_point(color="grey")+
  facet_wrap(~name,nrow=1,scales="free_y") +
  xlab("t") +
  ylab("Log expression") +
  theme_bw()

rownames(my.svg$f) <- rownames(Y.sub)
expr <- t(my.svg$f[c(gene.1,gene.2),]) |>
  as.data.frame() |>
  mutate(x=cso$t) |>
  pivot_longer(cols=-x)

plot2 <- expr |> ggplot(aes(x=x,y=exp(value))) +
  geom_line()+
  geom_hline(yintercept = 1, color="red", linetype="dashed")+
  facet_wrap(~name,nrow=1) +
  xlab("t") +
  ylab("Fold change from baseline") +
  theme_bw()

plot.big <- ggarrange(plot1,
          fit$plot,
          plot2,
          nrow=3,ncol=1)

rownames(my.svg$res) <- rownames(Y.sub)

ixs <- which(my.svg$res$p.val < 0.05/nrow(Y.sub))

peak.retained <- my.svg$res[ixs,]

peak.retained <- peak.retained |> as.data.frame() |>
  arrange(desc(peak))

peaks <- my.svg$res |> as.data.frame() |>
  arrange(desc(peak))

ixs <- c("Nlrc5", "Ddx58", "Dhx58", "Tnfsf10", "Ifih1")

expr <- t(Y.sub[ixs,]) |>
  apply(2, function(x) x/max(x)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))

plot.cs <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.15, alpha=0.6) +
  scale_color_gradient(low="grey90", high="darkblue")+
  facet_wrap(~name, ncol=1) +
  theme_bw()

ggsave(plot.cs, filename="../plots/top6_cs.png", width=5, height=20,units="in")


gene.1 <- which(rownames(Y.sub) == "Ifih1")
gene.2 <- which(rownames(Y.sub) == "Kitl")

#Plot these two
expr <- t(Y.sub[c(gene.1,gene.2),]) |>
  apply(2,function(x) x/max(x)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))
plot1 <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.25,alpha=0.75) +
  scale_color_gradient(low="grey90", high="darkblue")+
  facet_wrap(~name) +
  theme_bw()



gene.1 <- "P2ry12"
ixs <- which(t > 0.3 & t < 0.4)
expr <- t(Y.sub[gene.1,ixs]) |>
  apply(2,function(x) log(x+1)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x[ixs], y=meta.sub$y[ixs]) |>
  pivot_longer(cols=-c(x,y))

plot1 <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point(size=0.25,alpha=0.75) +
  scale_color_gradient(low="grey90", high="darkblue")+
  facet_wrap(~name) +
  theme_bw()


rownames(my.svg$f) <- rownames(Y.sub)

ixs <- which(apply(my.svg$f, 1, max) > 0.01)

F <- my.svg$f[ixs,]

my.pca <- prcomp(F)


