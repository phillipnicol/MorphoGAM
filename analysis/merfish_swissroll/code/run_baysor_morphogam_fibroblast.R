setwd(here::here("analysis/merfish_swissroll/code"))

### Curve1

load("../data/061923_D9_m2_Swiss_fibroblast.RData")

Y <- Matrix::readMM("../data/061923_D9_m2_Swiss_fibroblast_counts.mtx")

Y <- as.matrix(Y)


mgam <- MorphoGAM(Y=Y,curve.fit=fit,
                  design=y~s(t) + s(r))

save(mgam, file="../data/mgam_d9_061923_fibroblast_baysor.RData")


### Curve 2


load("../data/080823_D9_m5_Swiss_fibroblast.RData")

Y <- Matrix::readMM("../data/080823_D9_m5_Swiss_fibroblast_counts.mtx")

Y <- as.matrix(Y)


mgam <- MorphoGAM(Y=Y,curve.fit=fit,
                  design=y~s(t) + s(r))

save(mgam, file="../data/mgam_d9_m5_080823_fibroblast_baysor.RData")



### Curve 3

load("../data/080823_D9_m13_Swiss_fibroblast.RData")

Y <- Matrix::readMM("../data/080823_D9_m13_Swiss_fibroblast_counts.mtx.mtx")

Y <- as.matrix(Y)


mgam <- MorphoGAM(Y=Y,curve.fit=fit,
                  design=y~s(t) + s(r))

save(mgam, file="../data/mgam_d9_m13_080823_fibroblast_baysor.RData")

