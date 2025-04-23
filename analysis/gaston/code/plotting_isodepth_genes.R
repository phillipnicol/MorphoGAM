setwd(here::here("analysis/gaston/code"))


isodepth <- read.csv(file="../data/gaston_isodepth.csv")

isodepth <- as.vector(isodepth[,1])

plot(isodepth, Y.sub["Apob",])
