

### Curve1

load("../data/curve_d9_061923.RData")

Y <- Matrix::readMM("../data/d9_061923.mtx")

Y <- as.matrix(Y)


mgam <- MorphoGAM(Y=Y,curve.fit=fit,
                  design=y~s(t) + s(r))

save(mgam, file="../data/mgam_d9_061923_baysor.RData")


### Curve 2


load("../data/curve_d9_m5_080823.RData")

Y <- Matrix::readMM("../data/080823_d9_m5.mtx")

Y <- as.matrix(Y)


mgam <- MorphoGAM(Y=Y,curve.fit=fit,
                  design=y~s(t) + s(r))

save(mgam, file="../data/mgam_d9_m5_080823_baysor.RData")



### Curve 3

load("../data/curve_d9_m13_080823.RData")

Y <- Matrix::readMM("../data/080823_d9_m13.mtx")

Y <- as.matrix(Y)


mgam <- MorphoGAM(Y=Y,curve.fit=fit,
                  design=y~s(t) + s(r))

save(mgam, file="../data/mgam_d9_m13_080823_baysor.RData")

