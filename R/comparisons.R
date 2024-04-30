

library(ggplot2)
set.seed(1)
n <- 10^3

t <- runif(n=n, min=0, max=1)

r <- 0.5 + 0.5*t
xy <- data.frame(x=r*cos(10*t)+rnorm(n=n,sd=0.01),
                 y=r*sin(10*t) + rnorm(n=n,sd=0.01))

p <- xy |> ggplot(aes(x=x,y=y)) + geom_point(size=0.5) +
  theme_bw() +
  ggtitle("Spiral")

library(princurve)

fit.ps <- principal_curve(as.matrix(xy))
fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)

fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p1 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("Principal Curve")



fit.ps <- principal_curve(as.matrix(xy),df=50)
fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)

fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p2 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("Principal Curve")




