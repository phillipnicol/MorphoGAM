makeSwissRoll <- function(n=10^3,sd=0.03) {
  t <- runif(n=n, min=0, max=1)

  r <- 0.5 + 0.5*t
  xy <- data.frame(x=r*cos(10*t)+rnorm(n=n,sd=sd),
                   y=r*sin(10*t) + rnorm(n=n,sd=sd))
  return(as.matrix(xy))
}
