
spe <- STexampleData::Visium_humanDLPFC()

spe <- spe[,!is.na(spe$ground_truth)]
spe <- spe[,spe$ground_truth == "Layer4"]


radial.p <- rep(1, nrow(Y.sub))
angle.p <- rep(1, nrow(Y.sub))
angle <- rep(0,nrow(Y.sub))
radial <- rep(0,nrow(Y.sub))
angle.response <- angle
radial.response <- radial
l.o <- log(colSums(Y.sub))
fit <- CurveFinder(locus,knn=10,loop=TRUE)
fx.r <- Y.sub
fx.t <- Y.sub

Y <- counts(spe)
radial.p <- rep(1, nrow(Y))
angle.p <- rep(1, nrow(Y))
angle <- rep(0,nrow(Y))
radial <- rep(0,nrow(Y))
angle.response <- angle
radial.response <- radial
l.o <- log(colSums(Y))
fit <- CurveFinder(spatialCoords(spe))
fx.r <- Y
fx.t <- Y

for(i in 1:nrow(Y)) {
  print(i)
  fit1 <- mgcv::gam(Y[i,] ~ s(fit$xyt$t,bs="cr") + s(fit$xyt$r,bs="cr") + offset(l.o), family=nb(),
                    H = diag(19))

  se_beta <- sqrt(diag(fit1$rV %*% t(fit1$rV))[-1])

  #Shrink
  my.ash <- ashr::ash(fit1$coefficients[-1], se_beta)
  beta.shrink <- apply(ashr::get_post_sample(my.ash,1000),
                       2,
                       median)

  mat <- as.matrix(mgcv::predict.gam(fit1, type = "lpmatrix")[,-1])
  fx1 <- as.vector(mat[,1:9] %*% beta.shrink[1:9])
  fx2 <- as.vector(mat[,10:18] %*% beta.shrink[10:18])

  #fx.t[i,] <- fx1
  #fx.r[i,] <- fx2

  angle.response[i] <- max(exp(fit1$coefficients[1] + fx1)) - min(exp(fit1$coefficients[1] + fx1))
  radial.response[i] <- max(exp(fit1$coefficients[1] + fx2)) - min(exp(fit1$coefficients[1] + fx2))

  #pp <- predict(fit1, type="terms")

  angle[i] <- max(abs(fx1))
  angle.p[i] <- summary(fit1)$s.pv[1]
  radial[i] <- max(abs(fx2))
  radial.p[i] <- summary(fit1)$s.pv[2]
}


l.o <- log(colSums(Y))
for(i in 1:nrow(Y)) {
  if(sum(Y[i,]) <= 10) {
    next
  }
  print(i)
  fit1 <- mgcv::gam(Y[i,] ~ s(fit$xyt$t,bs="cr") + offset(l.o), family=nb(),
                    H = diag(10))

  se_beta <- sqrt(diag(fit1$rV %*% t(fit1$rV))[-1])

  #Shrink
  my.ash <- ashr::ash(fit1$coefficients[-1], se_beta)
  beta.shrink <- apply(ashr::get_post_sample(my.ash,1000),
                       2,
                       median)

  mat <- as.matrix(mgcv::predict.gam(fit1, type = "lpmatrix")[,-1])
  fx1 <- as.vector(mat[,1:9] %*% beta.shrink[1:9])
  #fx2 <- as.vector(mat[,10:18] %*% beta.shrink[10:18])

  #fx.t[i,] <- fx1
  #fx.r[i,] <- fx2

  angle.response[i] <- max(exp(fit1$coefficients[1] + fx1)) - min(exp(fit1$coefficients[1] + fx1))
  #radial.response[i] <- max(exp(fit1$coefficients[1] + fx2)) - min(exp(fit1$coefficients[1] + fx2))

  #pp <- predict(fit1, type="terms")

  angle[i] <- max(abs(fx1))
  angle.p[i] <- summary(fit1)$s.pv[1]
  #radial[i] <- max(abs(fx2))
  #radial.p[i] <- summary(fit1)$s.pv[2]
}





fit$curve.plot + ggtitle("Layer 4, DLPFC")







df <- data.frame(x=spatialCoords(spe.full)[,1],
                 y=spatialCoords(spe.full)[,2],
                 level = counts(spe.full)[17154,]) |>
  ggplot(aes(x=x,y=y,color=level)) +
  geom_point() +
  scale_color_gradient(low="grey", high="firebrick")



