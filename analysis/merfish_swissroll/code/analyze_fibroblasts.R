setwd(here::here("analysis/merfish_swissroll/code"))


load("../data/mgam_d9_061923_fibroblast_baysor.RData")

load("../data/061923_D9_m2_Swiss_fibroblast.RData")

Y <- Matrix::readMM("../data/061923_D9_m2_Swiss_fibroblast_counts.mtx")

gene_names <- readRDS("../data/gene_names.RDS")

library(MorphoGAM)


rownames(Y) <- gene_names
Y1 <- Y
rownames(mgam$results) <- gene_names
rownames(mgam$fxs.r) <- gene_names
rownames(mgam$fxs.t) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names

fit1 <- fit

# Find genes correlated with Vil1
vil1.1 <- mgam$fxs.t["Vil1",]

cos.sim1 <- apply(mgam$fxs.t, 1, function(x) {
  sum(x*vil1.1)/(sqrt(sum(x^2)*sum(vil1.1^2)))
})

mgam1 <- mgam


## Two


load("../data/mgam_d9_m5_080823_fibroblast_baysor.RData")

load("../data/080823_D9_m5_Swiss_fibroblast.RData")

Y <- Matrix::readMM("../data/080823_D9_m5_Swiss_fibroblast_counts.mtx")


rownames(Y) <- gene_names
Y2 <- Y
rownames(mgam$results) <- gene_names
rownames(mgam$fxs.r) <- gene_names
rownames(mgam$fxs.t) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names

fit2 <- fit

# Find genes correlated with Vil1
vil1.2 <- mgam$fxs.t["Vil1",]

cos.sim2 <- apply(mgam$fxs.t, 1, function(x) {
  sum(x*vil1.2)/(sqrt(sum(x^2)*sum(vil1.2^2)))
})

cos.sim2[order(abs(cos.sim2), decreasing = TRUE)[1:10]]

mgam2 <- mgam


## Third



load("../data/mgam_d9_m13_080823_fibroblast_baysor.RData")

load("../data/080823_D9_m13_Swiss_fibroblast.RData")

Y <- Matrix::readMM("../data/080823_D9_m13_Swiss_fibroblast_counts.mtx")




rownames(Y) <- gene_names
Y3 <- Y
rownames(mgam$results) <- gene_names
rownames(mgam$fxs.r) <- gene_names
rownames(mgam$fxs.t) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names
rownames(mgam$fpca.t$u) <- gene_names

fit3 <- fit

# Find genes correlated with Vil1
vil1.3 <- mgam$fxs.t["Vil1",]

cos.sim3 <- apply(mgam$fxs.t, 1, function(x) {
  sum(x*vil1.3)/(sqrt(sum(x^2)*sum(vil1.3^2)))
})

cos.sim3[order(abs(cos.sim3), decreasing = TRUE)[1:10]]

mgam3 <- mgam


df <- data.frame(roll1 = cos.sim1,
                 roll2 = cos.sim2,
                 roll3 = cos.sim3)


df$neg.markers <- pmin(df$roll1, df$roll2, df$roll3)
df$pos.markers <- pmax(df$roll1, df$roll2, df$roll3)

library(tidyverse)

neg.markers <- df |> arrange(desc(neg.markers)) |>
  select(roll1, roll2, roll3) |> head(n=6) |> rownames()

neg.markers <- neg.markers[-1]

pos.markers <- df |> arrange(pos.markers) |>
  select(roll1, roll2, roll3) |> head(n=5) |> rownames()

#pos.markers <- c("Syt11", "Fpr1", "Cr2", "Plek", "Mstn")

p1 <- t(mgam3$fxs.t[pos.markers,]) |> as.data.frame() |>
  mutate(t = fit3$xyt$t,) |>
  pivot_longer(cols=-t) |>
  ggplot(aes(x=t,y=value, color=name)) +
  geom_line() +
  theme_bw() +
  annotate("rect", xmin = 0.4, xmax = 0.7, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  labs(color="Gene")


p2 <- t(mgam2$fxs.t[pos.markers,]) |> as.data.frame() |>
  mutate(t = fit2$xyt$t) |>
  pivot_longer(cols=-t) |>
  mutate(min.fx = min(value), max.fx = max(value)) |>
  ggplot(aes(x=t,y=value, color=name)) +
  geom_line() +
  theme_bw() +
  annotate("rect", xmin = 0.25, xmax = 0.6, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  annotate("rect", xmin = 0.78, xmax = 0.8, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  annotate("rect", xmin = 0.88, xmax = 0.9, ymin = -Inf, ymax = Inf,
         fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  labs(color="Gene")



p3 <- t(mgam1$fxs.t[pos.markers,]) |> as.data.frame() |>
  mutate(t = fit1$xyt$t) |>
  pivot_longer(cols=-t) |>
  ggplot(aes(x=t,y=value, color=name)) +
  geom_line() +
  theme_bw() +
  annotate("rect", xmin = 0.95, xmax = 0.97, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  labs(color="Gene")

p.pos <- ggarrange(p1, p2, p3, nrow=3, ncol=1, common.legend=TRUE, legend="top")

#neg.markers <- c("Mafk", "Nlrp9b", "Fosl2", "Foxo3")




p1 <- t(mgam3$fxs.t[neg.markers,]) |> as.data.frame() |>
  mutate(t = fit3$xyt$t,) |>
  pivot_longer(cols=-t) |>
  ggplot(aes(x=t,y=value, color=name)) +
  geom_line() +
  theme_bw() +
  annotate("rect", xmin = 0.4, xmax = 0.7, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  labs(color="Gene")


p2 <- t(mgam2$fxs.t[neg.markers,]) |> as.data.frame() |>
  mutate(t = fit2$xyt$t) |>
  pivot_longer(cols=-t) |>
  mutate(min.fx = min(value), max.fx = max(value)) |>
  ggplot(aes(x=t,y=value, color=name)) +
  geom_line() +
  theme_bw() +
  annotate("rect", xmin = 0.25, xmax = 0.6, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  annotate("rect", xmin = 0.78, xmax = 0.8, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  annotate("rect", xmin = 0.88, xmax = 0.9, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  labs(color="Gene")



p3 <- t(mgam1$fxs.t[neg.markers,]) |> as.data.frame() |>
  mutate(t = fit1$xyt$t) |>
  pivot_longer(cols=-t) |>
  ggplot(aes(x=t,y=value, color=name)) +
  geom_line() +
  theme_bw() +
  annotate("rect", xmin = 0.95, xmax = 0.97, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  labs(color="Gene")


p.neg <- ggarrange(p1, p2, p3, nrow=3, ncol=1, common.legend=TRUE, legend="top")

ggsave(p.neg, filename="../plots/ulcerated_neg.png",
       width=7.32, height=8.36, units="in")

ggsave(p.pos, filename="../plots/ulcerated_pos.png",
       width=7.32, height=8.36, units="in")

