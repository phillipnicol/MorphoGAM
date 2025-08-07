setwd(here::here("analysis/sim_curve/code"))

library(ggplot2)

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

n <- 1000

r_inner <- 1
r_outer <- 1.25
xyr <- generate_annulus(n,r_inner=1,r_outer=1.25)

library(MorphoGAM)

fit <- MorphoGAM::CurveFinder(xy=cbind(xyr$x, xyr$y) |> as.matrix(), loop=TRUE,k=10)

#Rotate the curve to start near (1,0)

fit$xyt$t <- fit$xyt$t * 2*pi
ix <- which.min(fit$xyt$t)
t.min <- min(fit$xyt$t)

#xyr$theta[ix]

Rotation <- (2*pi - xyr$theta[ix])

#t.new <- ifelse(fit$xyt$t + Rotation < 2*pi, fit$xyt$t + Rotation, fit$xyt$t + Rotation - 2*pi)
theta.new <- ifelse(xyr$theta - xyr$theta[ix] < 0, xyr$theta - xyr$theta[ix] + 2*pi, xyr$theta - xyr$theta[ix])

xyr$theta.new <- theta.new
p.true.r <- ggplot(data=xyr, aes(x=x,y=y,color=(r-r_inner)/(r_outer-r_inner))) + geom_point() + theme_bw() +
  scale_color_gradient(low="blue", high="red") + guides(color="none")
p.true.t <- ggplot(data=xyr, aes(x=x,y=y,color=theta)) + geom_point() + theme_bw() + scale_color_gradient(low="blue", high="red") +
  guides(color="none")

# Plot about t (theta)
df.t <- data.frame(x = xyr$theta.new, y = fit$xyt$t)
r2.t <- summary(lm(y ~ x, data = df.t))$r.squared

p.2 <- ggplot(data = df.t, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = Inf, y = -Inf, label = paste0("R² = ", round(r2.t, 3)),
           hjust = 1.1, vjust = -1.1, size = 4) +
  theme_bw() + xlab("True t") + ylab("Estimated t")

# Plot about r (radius)
df.r <- data.frame(x = (xyr$r - r_inner)/(r_outer - r_inner), y = fit$xyt$r)
r2.r <- summary(lm(y ~ x, data = df.r))$r.squared

p.3 <- ggplot(data = df.r, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = Inf, y = -Inf, label = paste0("R² = ", round(r2.r, 3)),
           hjust = 1.1, vjust = -1.1, size = 4) +
  theme_bw() + xlab("True r") + ylab("Estimated r")

library(ggpubr)

p1 <- ggarrange(p.true.t, p.true.r, p.2, p.3,
          nrow = 1, labels=c("a","b","c","d"))


p1 <- annotate_figure(p1, top=text_grob("Annulus width = 0.25", size=14))




annulus_performance <- function(k, width, n=1000) {
  r_inner <- 1
  r_outer <- 1 + width
  xyr <- generate_annulus(n,r_inner=r_inner,r_outer=r_outer)

  library(MorphoGAM)

  fit <- MorphoGAM::CurveFinder(xy=cbind(xyr$x, xyr$y) |> as.matrix(), loop=TRUE,k=k)

  #Rotate the curve to start near (1,0)



  #Rotate the curve to start near (1,0)

  fit$xyt$t <- fit$xyt$t * 2*pi
  ix <- which.min(fit$xyt$t)
  t.min <- min(fit$xyt$t)

  #xyr$theta[ix]

  Rotation <- (2*pi - xyr$theta[ix])

  #t.new <- ifelse(fit$xyt$t + Rotation < 2*pi, fit$xyt$t + Rotation, fit$xyt$t + Rotation - 2*pi)
  theta.new <- ifelse(xyr$theta - xyr$theta[ix] < 0, xyr$theta - xyr$theta[ix] + 2*pi, xyr$theta - xyr$theta[ix])

  xyr$theta.new <- theta.new

  t_r2 <- cor(xyr$theta.new, fit$xyt$t)^2

  r_r2 <- cor(xyr$r, fit$xyt$r)^2

  return(c(t_r2, r_r2))

}


set.seed(44)

param.test <- expand.grid(k=c(10,15,20), width=c(0.1,0.25,0.5,1,2.5,5))
param.test$t_r2 <- 0
param.test$r_r2 <- 0

reps <- 10
for(i in 1:nrow(param.test)) {
  print(i)
  for(j in 1:10) {
    v <- annulus_performance(k=param.test$k[i], width=param.test$width[i], n=1000)
    param.test$t_r2[i] <- param.test$t_r2[i] + v[1]/reps
    param.test$r_r2[i] <- param.test$r_r2[i] + v[2]/reps
  }
}

param.test$k <- as.factor(param.test$k)
p.t <- ggplot(data=param.test,aes(x=width, y=t_r2, color=k, group=k)) +
            geom_point() + geom_line() + theme_bw() + scale_x_log10() +
  xlab("Annulus width") + ylab("R2 between true vs estimated t")

p.r <- ggplot(data=param.test,aes(x=width, y=r_r2, color=k, group=k)) +
  geom_point() + geom_line() + theme_bw() + scale_x_log10() +
  xlab("Annulus width") + ylab("R2 between true vs estimated r")


p.full <- ggarrange(p1, ggarrange(p.t,p.r, nrow=1,labels=c("e","f")),
                    heights=c(1.5,1),nrow=2)

ggsave(p.full, filename="../plots/simulated_annulus_full.png",
       width=12.1, height=7.34, units="in")


