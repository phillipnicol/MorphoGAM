

#indic <- ifelse(fit2$xyt$t > 0.25 & fit2$xyt$t < 0.6, 1,
#                ifelse(fit2$xyt$t > 0.78 & fit2$xyt$t < 0.8, 1,
#                       ifelse(fit2$xyt$t < 0.88 & fit2$xyt$t < 0.9, 1, 0)))


indic <- rep(0, length(fit2$xyt$t))

indic[fit2$xyt$t > 0.25 & fit2$xyt$t < 0.6] <- 1
indic[fit2$xyt$t > 0.78 & fit2$xyt$t < 0.8] <- 1
indic[fit2$xyt$t > 0.88 & fit2$xyt$t < 0.9] <- 1

fit <- mgcv::gam(formula = Y2[1,] ~ indic + s(fit2$xyt$r, bs="cr"),
                 family = mgcv::nb())

l.o <- log(colSums(as.matrix(Y2)))
res <- apply(Y2, 1, function(x) {
  fit <- mgcv::gam(formula = x ~ indic + s(fit2$xyt$r, bs="cr"),
                   family = mgcv::nb(),
                   offset = l.o)

  return(c(fit$coefficients[2], summary(fit)$pTerms.table[1,3]))
})





