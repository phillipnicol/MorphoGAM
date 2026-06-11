
set.seed(1)
n <- 1000
#Generate points on a square 
x <- runif(n, 0, 1)
y <- runif(n, 0, 1)

xy <- cbind(x, y)

fit <- CurveFinder(xy)