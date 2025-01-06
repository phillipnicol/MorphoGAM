
## Activate the script
#renv::activate()


here("analysis","merfish_swissroll","code") |> setwd()

#setwd("analysis/merfish_swissroll/code") #Assume MorphoGAM is current directory

Y <- read.csv("../../data/080823_D9_m5_Swiss.csv")

#Subample to hand-draw curve
set.seed(1)
Y.sub <- Y[sample(1:nrow(Y), 10^4,replace=FALSE),]

xy <- Y.sub[,c("abs_position_1", "abs_position_2")]

library(MorphoGAM)

xy <- as.matrix(xy)

fit <- CurveFinderInteractive(xy)

save(fit, file="../data/080823_D9_m5_Swiss_curve.RData")

xy.all <- Y[,c("abs_position_1", "abs_position_2")]
xy.all <- as.matrix(xy.all)

fit$xyt <- fit$xyt |> arrange(t)

t.all <- princurve::project_to_curve(x=xy.all,
                                     s=cbind(fit$xyt$f1,
                                             fit$xyt$f2))

t.all2 <- as.numeric(t.all$lambda/max(t.all$lambda))

saveRDS(t.all2, "../data/080823_D9_m5_Swiss_curve_t.RDS")

# Load the curve

load("../data/080823_D9_m5_Swiss_curve.RData")

my.t <- readRDS("../data/080823_D9_m5_Swiss_curve_t.RDS")

Y$t <- round(my.t, 4)

Y <- Y[,c("gene_name", "t")]

result_matrix <- Y %>%
  mutate(t = as.factor(t)) %>%
  group_by(gene_name, t) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = t, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "gene_name")


result_matrix <- as.matrix(result_matrix)

#t.use <- round(my.t,4)
l.o <- log(colSums(result_matrix))
fxs.t <- matrix(0,nrow=nrow(result_matrix), ncol=ncol(result_matrix))
for(i in 1:nrow(result_matrix)) {
  print(i)
  gam.fit <- mgcv::gam(result_matrix[i,] ~ s(colnames(result_matrix) |> as.numeric(), bs="cr"),
                       family=nb(), offset=l.o)

  basis.functions <- mgcv::predict.gam(gam.fit, type="lpmatrix")[,-1]
  fxs.t[i,] <- basis.functions %*% coefficients(gam.fit)[-1]
}


mgam <- list()
mgam$fxs.t <- fxs.t
rownames(mgam$fxs.t) <- rownames(result_matrix)
mgam$fpca <- irlba::irlba(fxs.t)

num_genes <- 5

top5 <- order(abs(mgam$fpca$u[,L]), decreasing = TRUE)[1:num_genes]

if(sum(mgam$fxs.t[top5,]%*%mgam$fpca$v[,L]) < 0) {
  mgam$fpca$v[,L] <- -1*mgam$fpca$v[,L]
}

mat <- matrix(0, nrow=ncol(mgam$fxs.t), ncol=num_genes + 1)
mat[,1] <- max(mgam$fxs.t[top5,])*(mgam$fpca$v[,L]/max(mgam$fpca$v[,L])) #Max 1
colnames(mat) <- rep("Eigenfunction", num_genes+1)
for(i in 1:num_genes) {
  mat[,i+1] <- mgam$fxs.t[top5[i],]
  colnames(mat)[i+1] <- rownames(mgam$fxs.t)[top5[i]]
}
#mat[,1] <- mgam$fpca$v[,L]*max(mgam$fxs.t[top5,])

df <- as.data.frame(mat)
df$t <- colnames(result_matrix) |> as.numeric()
df <- df |> pivot_longer(cols=-t)

p <- ggplot(data=df,aes(x=t,y=value, group=name)) +
  geom_line()

# Generate a bright palette that avoids yellow
#num_genes <- length(unique(df$name))
base_colors <- brewer.pal(9, "Set1")[-6]  # Exclude the 6th color (yellow)

# Extend the palette dynamically if more colors are needed
if (num_genes <= length(base_colors)) {
  custom_colors <- base_colors[1:num_genes]
} else {
  custom_colors <- colorRampPalette(base_colors)(num_genes)
}

# Ensure "Eigenfn" is black
names(custom_colors) <- colnames(mat)[-1]
custom_colors["Eigenfunction"] <- "black"

label_data <- df %>%
  group_by(name) %>%
  filter(value == max(value))  # Select rows where the value is maximum for each name


# Plot
p <- ggplot(data = df, aes(x = t, y = value,
                           group = name,
                           color = name)) +
  geom_line(aes(size = ifelse(name == "Eigenfunction", 2, 0.75)),
            show.legend = FALSE) +  # Conditional size
  scale_color_manual(values = custom_colors) +  # Custom color palette
  scale_size_identity() +  # Use size as is, without scaling
  theme_bw() +
  geom_text_repel(data = label_data,
                  aes(label = name),
                  size = 4,
                  nudge_x = 0.1,
                  hjust = 0,
                  show.legend = FALSE) +
  ylab("Log FC from baseline") +
  xlab("t") +
  geom_hline(yintercept = 0, linetype="dashed",
             color="blue")





## Bin based on x and y
Y.sub <- Y.sub[,c("gene_name", "abs_position_1", "abs_position_2")]
Y.sub <- as.data.frame(Y.sub)
colnames(Y.sub) <- c("gene","x","y")

Y.sub$x <- (Y.sub$x - min(Y.sub$x))/(max(Y.sub$x) - min(Y.sub$x))
Y.sub$y <- (Y.sub$y - min(Y.sub$y))/(max(Y.sub$y) - min(Y.sub$y))

# Discretize
Y.sub$x <- round(Y.sub$x,digits=2)
Y.sub$y <- round(Y.sub$y, digits=2)

Y.sub$xy <- paste(Y.sub$x, Y.sub$y, sep = ",")
gene_matrix <- table(Y.sub$gene, Y.sub$xy)
gene_matrix <- gene_matrix[,colSums(gene_matrix) > 210]

# Extract x and y coordinates from the column names of gene_matrix
coordinates <- strsplit(colnames(gene_matrix), ",")
coordinates_matrix <- do.call(rbind, coordinates)
colnames(coordinates_matrix) <- c("x", "y")
#coordinates_matrix <- as.data.frame(coordinates_matrix)

# Convert x and y to numeric for clarity (optional)
mode(coordinates_matrix) <- "numeric"

#xy <- cbind(Y.sub$x, Y.sub$y)

xy <- coordinates_matrix
fit <- CurveFinderInteractive(xy)

mgam <- MorphoGAM(gene_matrix, curve.fit=fit,design=y~s(t)+s(r))




