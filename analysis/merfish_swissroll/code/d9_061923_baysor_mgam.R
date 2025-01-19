

#library(MorphoGAM)

#xy <- matrix(cbind(ad$obs$x, ad$obs$y),ncol=2)


#fit <- CurveFinderInteractive(xy)
#
#save(fit, file="../data/curve_d9_061923_baysor.RData")

Y <- ad$raw$X |> t()

rownames(Y) <- ad$var_names


mgam <- MorphoGAM(Y=Y,curve.fit=fit,
                  design=y~s(t) + s(r))
