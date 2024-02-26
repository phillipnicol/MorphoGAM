

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

gene.1 <- which(rownames(Y.sub) == "Kitl")
gene.2 <- which(rownames(Y.sub) == "Tnfsf10")

#Plot these two
expr <- t(Y.sub[c(gene.1,gene.2),]) |>
  apply(2,function(x) x/max(x)) |>
  as.data.frame() |>
  mutate(x=meta.sub$x, y=meta.sub$y) |>
  pivot_longer(cols=-c(x,y))

plot <- expr |> ggplot(aes(x=x,y=y,color=value)) +
  geom_point() +
  scale_color_gradient(low="grey90", high="darkblue")+
  facet_wrap(~name) +
  theme_bw()

#Test expression for these genes
library(SPARK)
locus <- as.matrix(meta.sub[,c("x","y")])
res <- SPARK::sparkx(count_in = Y.sub,
                     locus_in = locus)

p.vals <- res$res_mtest |> as.data.frame() |>
  arrange(adjustedPval)



