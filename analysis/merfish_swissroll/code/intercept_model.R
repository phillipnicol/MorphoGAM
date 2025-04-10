
### Roll 1

indic <- rep(0, length(fit1$xyt$t))

indic[fit1$xyt$t > 0.95 & fit1$xyt$t < 0.97] <- 1

l.o <- log(colSums(as.matrix(Y1)))
res1 <- apply(Y1, 1, function(x) {
  fit <- mgcv::gam(formula = x ~ indic + s(fit1$xyt$r, bs="cr"),
                   family = mgcv::nb(),
                   offset = l.o)

  return(c(fit$coefficients[2], summary(fit)$pTerms.table[1,3]))
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

  return(c(fit$coefficients[2], summary(fit)$pTerms.table[1,3]))
})



### Roll 3

indic <- rep(0, length(fit3$xyt$t))

indic[fit3$xyt$t > 0.4 & fit3$xyt$t < 0.7] <- 1

l.o <- log(colSums(as.matrix(Y3)))
res3 <- apply(Y3, 1, function(x) {
  fit <- mgcv::gam(formula = x ~ indic + s(fit3$xyt$r, bs="cr"),
                   family = mgcv::nb(),
                   offset = l.o)

  return(c(fit$coefficients[2], summary(fit)$pTerms.table[1,3]))
})




## Analysis

res1 <- t(res1)
res2 <- t(res2)
res3 <- t(res3)




library(tidyverse)
df <- data.frame(gene = rownames(res1),
                 beta1 = res1[,1],
                 beta2 = res2[,1],
                 beta3 = res3[,1],
                 pval = nrow(res1)*pmax(res1[,2], res2[,2], res3[,2]))

#df <- df |> mutate(stat = -2*stat) |>
#  mutate(pval = pchisq(stat, df = 6, lower.tail=FALSE)) |>
#  mutate(pval = pmin(pval*nrow(res1), 1))

#df <- df |> mutate(pval = round(pmin(pval, 1), digits=3))

genes <- df |> mutate(max.min = pmin(beta1, beta2, beta3)) |>
  arrange(desc(max.min)) |> head(n=20) |> select(gene)

#df |> arrange(pval) |> head(n=20)




genes <- df |> mutate(min.max = pmax(beta1, beta2, beta3)) |>
  arrange(min.max) |> head(n=20) |> select(gene)

