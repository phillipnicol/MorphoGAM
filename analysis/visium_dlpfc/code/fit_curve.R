
setwd(here::here("analysis", "visium_dlpfc", "code"))

library(STexampleData)


spe <- Visium_humanDLPFC()

for(i in 2:6) {
  print(i)
  layer <- paste0("Layer", i)
  ixs <- which(spe$ground_truth == layer) #subset to Layer

  xy <- spatialCoords(spe)[ixs,]

  fit <- CurveFinder(xy, knn=10)

  mgam <- MorphoGAM(Y = counts(spe)[,ixs] |> as.matrix(),
                    curve.fit = fit,
                    design = y ~ s(t, bs="cr", k=10) + s(r, bs="cr", k=10))

  my.t <- mgam$results |> arrange(desc(peak.t))

  plotGAMestimates(Y, genes=rownames(my.t)[1:5], curve_fit=fit, mgam_object = mgam,nrow=1)


}
