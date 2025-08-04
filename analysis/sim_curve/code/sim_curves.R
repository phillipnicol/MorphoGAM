setwd(here::here("analysis/sim_curve/code"))

set.seed(123)  # for reproducibility

generate_annulus <- function(n=1000, r_inner, r_outer) {
  theta <- runif(n, 0, 2 * pi)         # angle uniformly from 0 to 2pi
  r_unif <- runif(n, 0, 1)             # uniform in [0,1] for radius transformation

  # To get uniform area density in the annulus, use square root transformation
  r <- sqrt(r_unif * (r_outer^2 - r_inner^2) + r_inner^2)

  # Convert polar to Cartesian coordinates
  x <- r * cos(theta)
  y <- r * sin(theta)

  data.frame(x = x, y = y, r=r, theta=theta)
}

#Case 1

r_inner <- 1
r_outer <- 1.25
xyr <- generate_annulus(n,r_inner=1,r_outer=1.25)

library(MorphoGAM)

fit <- MorphoGAM::CurveFinder(xy=cbind(xyr$x, xyr$y) |> as.matrix(), loop=TRUE)

#Rotate the curve to start near (1,0)

fit$xyt$t <- fit$xyt$t * 2*pi
ix <- which.min(fit$xyt$t)
t.min <- min(fit$xyt$t)

#xyr$theta[ix]

Rotation <- (2*pi - xyr$theta[ix])

t.new <- ifelse(fit$xyt$t + Rotation < 2*pi, fit$xyt$t + Rotation, fit$xyt$t + Rotation - 2*pi)

p.1 <- ggplot(data=xyr, aes(x=x,y=y,color=(r-r_inner)/(r_outer-r_inner))) + geom_point() + theme_bw()

#Plot about t

df <- data.frame(x=xyr$theta, y=t.new)
p.2 <- ggplot(data=df,aes(x=x,y=y)) + geom_point()

#Plot about r

df <- data.frame(x=(xyr$r - r_inner)/(r_outer - r_inner), y=fit$xyt$r)

p.3 <- ggplot(data=df,aes(x=x,y=y)) + geom_point()
