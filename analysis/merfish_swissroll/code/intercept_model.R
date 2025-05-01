
### Roll 1

indic <- rep(0, length(fit1$xyt$t))

indic[fit1$xyt$t > 0.95 & fit1$xyt$t < 0.97] <- 1

l.o <- log(colSums(as.matrix(Y1)))
res1 <- apply(Y1, 1, function(x) {
  fit <- mgcv::gam(formula = x ~ indic + s(fit1$xyt$r, bs="cr"),
                   family = mgcv::nb(),
                   offset = l.o)

  return(c(fit$coefficients[2],
           fit$coefficients[1],
           summary(fit)$pTerms.table[1,3]))
})






### Roll 2

#indic <- ifelse(fit2$xyt$t > 0.25 & fit2$xyt$t < 0.6, 1,
#                ifelse(fit2$xyt$t > 0.78 & fit2$xyt$t < 0.8, 1,
#                       ifelse(fit2$xyt$t < 0.88 & fit2$xyt$t < 0.9, 1, 0)))


indic <- rep(0, length(fit2$xyt$t))

indic[fit2$xyt$t > 0.25 & fit2$xyt$t < 0.6] <- 1
indic[fit2$xyt$t > 0.78 & fit2$xyt$t < 0.8] <- 1
indic[fit2$xyt$t > 0.88 & fit2$xyt$t < 0.9] <- 1

l.o <- log(colSums(as.matrix(Y2)))
res2 <- apply(Y2, 1, function(x) {
  fit <- mgcv::gam(formula = x ~ indic + s(fit2$xyt$r, bs="cr"),
                   family = mgcv::nb(),
                   offset = l.o)

  return(c(fit$coefficients[2],
           fit$coefficients[1],
           summary(fit)$pTerms.table[1,3]))
})



### Roll 3

indic <- rep(0, length(fit3$xyt$t))

indic[fit3$xyt$t > 0.4 & fit3$xyt$t < 0.7] <- 1

l.o <- log(colSums(as.matrix(Y3)))
res3 <- apply(Y3, 1, function(x) {
  fit <- mgcv::gam(formula = x ~ indic + s(fit3$xyt$r, bs="cr"),
                   family = mgcv::nb(),
                   offset = l.o)

  return(c(fit$coefficients[2],
           fit$coefficients[1],
           summary(fit)$pTerms.table[1,3]))
})




## Analysis

res1 <- t(res1)
res2 <- t(res2)
res3 <- t(res3)

saveRDS(res1, file = "../data/roll1_intercept_model.RDS")
saveRDS(res2, file = "../data/roll2_intercept_model.RDS")
saveRDS(res3, file = "../data/roll3_intercept_model.RDS")


res1 <- readRDS("../data/roll1_intercept_model.RDS")
res2 <- readRDS("../data/roll2_intercept_model.RDS")
res3 <- readRDS("../data/roll3_intercept_model.RDS")


library(tidyverse)
beta.combine <- data.frame(gene = rownames(res1),
                 beta1 = res1[,1],
                 beta2 = res2[,1],
                 beta3 = res3[,1],
                 pval = 1 - (1-pmax(res1[,3], res2[,3], res3[,3])^3))

beta.combine$pval <- p.adjust(beta.combine$pval, method = "BH")

#df <- df |> mutate(stat = -2*stat) |>
#  mutate(pval = pchisq(stat, df = 6, lower.tail=FALSE)) |>
#  mutate(pval = pmin(pval*nrow(res1), 1))

#df <- df |> mutate(pval = round(pmin(pval, 1), digits=3))

#up.genes <- df |> mutate(max.min = pmin(beta1, beta2, beta3)) |>
#  arrange(desc(max.min)) |> head(n=10) |> select(gene) |> unlist()

#df |> arrange(pval) |> head(n=20)

up.genes <- beta.combine |> mutate(max.min = pmin(beta1, beta2, beta3)) |>
  arrange(desc(max.min))


#bottom.genes <- df |> mutate(min.max = pmax(beta1, beta2, beta3)) |>
#  arrange(min.max) |> head(n=20) |> select(gene)



bottom.genes <- df |> mutate(min.max = beta1 + beta2 + beta3) |>
  arrange(min.max) |> head(n=20) |> select(gene)

fxs.t <- Y1 |> as.matrix()
### Roll 1

indic <- rep(0, length(fit1$xyt$t))

indic[fit1$xyt$t > 0.95 & fit1$xyt$t < 0.97] <- 1

fxs.t[,indic == 0] <-0
fxs.t[,indic == 1] <- res1[,1]


df <- t(fxs.t[up.genes,]) |>
  as.data.frame() |>
  mutate(t = fit1$xyt$t) |>
  pivot_longer(-t)

# Arrange and identify horizontal segments
segments <- df %>%
  arrange(name, t) %>%
  group_by(name) %>%
  mutate(next_t = lead(t),
         next_value = lead(value)) %>%
  filter(!is.na(next_t), value == next_value) %>%
  ungroup()




# Plot horizontal segments only
p1 <- ggplot(segments, aes(x = t, xend = next_t, y = value, yend = value, color = name)) +
  theme_bw() +
  annotate("rect", xmin = 0.95, xmax = 0.97, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  geom_segment(size = 1, alpha=1) +
  labs(color = "Gene") +
  xlim(c(0.9, 1))







fxs.t <- Y2 |> as.matrix()
### Roll 1

indic <- rep(0, length(fit2$xyt$t))

indic[fit2$xyt$t > 0.25 & fit2$xyt$t < 0.6] <- 1
indic[fit2$xyt$t > 0.78 & fit2$xyt$t < 0.8] <- 1
indic[fit2$xyt$t > 0.88 & fit2$xyt$t < 0.9] <- 1

fxs.t[,indic == 0] <- 0
fxs.t[,indic == 1] <- res2[,1]


df <- t(fxs.t[up.genes,]) |>
  as.data.frame() |>
  mutate(t = fit2$xyt$t) |>
  pivot_longer(-t)

# Arrange and identify horizontal segments
segments <- df %>%
  arrange(name, t) %>%
  group_by(name) %>%
  mutate(next_t = lead(t),
         next_value = lead(value)) %>%
  filter(!is.na(next_t), value == next_value) %>%
  ungroup()


# Plot horizontal segments only
p2 <- ggplot(segments, aes(x = t, xend = next_t, y = value, yend = value, color = name)) +
  theme_bw() +
  annotate("rect", xmin = 0.25, xmax = 0.6, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  annotate("rect", xmin = 0.78, xmax = 0.8, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  annotate("rect", xmin = 0.88, xmax = 0.9, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  geom_segment(size = 1, alpha=1) +
  labs(color = "Gene")








fxs.t <- Y3 |> as.matrix()
### Roll 3

indic <- rep(0, length(fit3$xyt$t))

indic[fit3$xyt$t > 0.4 & fit3$xyt$t < 0.7] <- 1

fxs.t[,indic == 0] <-0
fxs.t[,indic == 1] <- res3[,1]


df <- t(fxs.t[up.genes,]) |>
  as.data.frame() |>
  mutate(t = fit3$xyt$t) |>
  pivot_longer(-t)

# Arrange and identify horizontal segments
segments <- df %>%
  arrange(name, t) %>%
  group_by(name) %>%
  mutate(next_t = lead(t),
         next_value = lead(value)) %>%
  filter(!is.na(next_t), value == next_value) %>%
  ungroup()




# Plot horizontal segments only
p3 <- ggplot(segments, aes(x = t, xend = next_t, y = value, yend = value, color = name)) +
  theme_bw() +
  annotate("rect", xmin = 0.4, xmax = 0.7, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  geom_segment(size = 1, alpha=1) +
  labs(color = "Gene")

library(ggpubr)

p.pos <- ggarrange(p1, p2, p3, nrow=3, ncol=1, common.legend=TRUE, legend="top")

ggsave(p.pos, filename="../plots/ulcerated_intercept_pos.png",
       width=7.32, height=8.36, units="in")







### Counts

#gene.names <- c("Mmp3", "Selp", "Mmp10", "Col12a1", "Il11")
#gene.names <- c("Selp", "Ghsr", "Blank-32", "Nlrp9a", "Plek")
gene.names <- c("Igfbp6", "Wnt2b", "Ackr4", "Dkk3")


Y1.dn <- sweep(as.matrix(Y1), MARGIN = 2, FUN = "/", STATS = colSums(Y1))
fxs.t <- Y1
indic <- rep(0, length(fit1$xyt$t))

indic[fit1$xyt$t > 0.95 & fit1$xyt$t < 0.97] <- 1

fxs.t[,indic == 0] <- exp(res1[,2])
fxs.t[,indic == 1] <- exp(res1[,2] + res1[,1])

df <- t(Y1.dn[gene.names,]) |>
  as.data.frame() |>
  mutate(t = fit1$xyt$t) |>
  pivot_longer(-t)


df.fit <- t(fxs.t[gene.names,]) |>
  as.data.frame() |>
  mutate(t = fit1$xyt$t) |>
  pivot_longer(-t)

df$fitted <- df.fit$value

# Arrange and identify horizontal segments
segments <- df %>%
  arrange(name, t) %>%
  group_by(name) %>%
  mutate(next_t = lead(t),
         next_value = lead(value)) %>%
  filter(!is.na(next_t), value == next_value) %>%
  ungroup()




# Plot horizontal segments only
p1 <- ggplot(df, aes(x = t, y = value)) +
  theme_bw() +
  geom_point(size=0.25) +
  annotate("rect", xmin = 0.95, xmax = 0.97, ymin = 0, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  facet_wrap(~name, ncol= 1) +
  geom_line(aes(x = t, y = fitted, color = name)) +
  labs(color = "Gene") +
  scale_y_sqrt()



### Roll 2


Y2 <- as.matrix(Y2)
Y2.dn <- sweep(as.matrix(Y2), MARGIN = 2, FUN = "/", STATS = colSums(Y2))
fxs.t <- Y2 |> as.matrix()
indic <- rep(0, length(fit2$xyt$t))


indic[fit2$xyt$t > 0.25 & fit2$xyt$t < 0.6] <- 1
indic[fit2$xyt$t > 0.78 & fit2$xyt$t < 0.8] <- 1
indic[fit2$xyt$t > 0.88 & fit2$xyt$t < 0.9] <- 1

fxs.t[,indic == 0] <- exp(res2[,2])
fxs.t[,indic == 1] <- exp(res2[,2] + res1[,1])

df <- t(Y2.dn[gene.names,]) |>
  as.data.frame() |>
  mutate(t = fit2$xyt$t) |>
  pivot_longer(-t)

df.fit <- t(fxs.t[gene.names,]) |>
  as.data.frame() |>
  mutate(t = fit2$xyt$t) |>
  pivot_longer(-t)

df$fitted <- df.fit$value

# Arrange and identify horizontal segments
segments <- df %>%
  arrange(name, t) %>%
  group_by(name) %>%
  mutate(next_t = lead(t),
         next_value = lead(value)) %>%
  filter(!is.na(next_t), value == next_value) %>%
  ungroup()




# Plot horizontal segments only
p2 <- ggplot(df, aes(x = t, y = value)) +
  theme_bw() +
  geom_point(size=0.25) +
  annotate("rect", xmin = 0.25, xmax = 0.6, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  annotate("rect", xmin = 0.78, xmax = 0.8, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  annotate("rect", xmin = 0.88, xmax = 0.9, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  geom_line(aes(x = t, y = fitted, color = name)) +
  facet_wrap(~name, ncol = 1) +
  labs(color = "Gene")



### Roll 3

indic <- rep(0, length(fit3$xyt$t))

indic[fit3$xyt$t > 0.4 & fit3$xyt$t < 0.7] <- 1


Y3 <- as.matrix(Y3)
Y3.dn <- sweep(as.matrix(Y3), MARGIN = 2, FUN = "/", STATS = colSums(Y3))
fxs.t <- Y3 |> as.matrix()

fxs.t[,indic == 0] <- exp(res3[,2])
fxs.t[,indic == 1] <- exp(res3[,2] + res3[,1])

df <- t(Y3.dn[gene.names,]) |>
  as.data.frame() |>
  mutate(t = fit3$xyt$t) |>
  pivot_longer(-t)

df.fit <- t(fxs.t[gene.names,]) |>
  as.data.frame() |>
  mutate(t = fit3$xyt$t) |>
  pivot_longer(-t)


df$fitted <- df.fit$value

# Plot horizontal segments only
p3 <- ggplot(df, aes(x = t, y = value)) +
  geom_point(size=0.25) +
  theme_bw() +
  annotate("rect", xmin = 0.4, xmax = 0.7, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5) +
  ylab("Log FC from baseline") +
  geom_line(aes(x = t, y = fitted, color = name)) +
  facet_wrap(~name, ncol = 1) +
  labs(color = "Gene")



