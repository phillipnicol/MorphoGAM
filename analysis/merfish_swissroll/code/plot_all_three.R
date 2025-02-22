

load("../data/curve_d9_m13_080823.RData")

p1 <- fit$curve.plot + guides(color="none")

load("../data/mgam_d9_m13_080823.RData")

mgam1 <- mgam


load("../data/curve_d9_080823.RData")

p2 <- fit$curve.plot + guides(color="none")

load("../data/mgam_d9_080823.RData")

mgam2 <- mgam

load("../data/curve_d9_061923.RData")

p3 <- fit$curve.plot + guides(color="none")

load("../data/mgam_d9_061923.RData")

mgam3 <- mgam

ggsave(ggarrange(p1,p2,p3,nrow=1),
       filename="all_swissroll.png")

