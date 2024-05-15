

library(ggplot2)
set.seed(1)
n <- 10^3

t <- runif(n=n, min=0, max=1)

r <- 0.5 + 0.5*t
xy <- data.frame(x=r*cos(10*t)+rnorm(n=n,sd=0.03),
                 y=r*sin(10*t) + rnorm(n=n,sd=0.03))

p <- xy |> ggplot(aes(x=x,y=y)) + geom_point(size=0.5) +
  theme_bw() +
  ggtitle("Spiral") + coord_fixed()

fit <- CurveFinder(xy |> as.matrix(), knn=10)

p.cp <- fit$curve.plot + guides(color="none") + ggtitle("Fitted curve")

p.p2 <- fit$coordinate.plot + guides(color="none")

p.p3 <- fit$residuals.plot + guides(color="none")



library(princurve)

fit.ps <- principal_curve(as.matrix(xy))
fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)

fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p.pc.df <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("Principal Curve") + guides(color="none")


fit.ps <- principal_curve(as.matrix(xy),df=50)
fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)

p.pc.df25 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("Principal Curve (df=35)") +guides(color="none")

fit.ps <- principal_curve(as.matrix(xy),df=100)
fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)

fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p.pc.df50 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("Principal Curve (df=50)") +guides(color="none")



knn.seq <- c(5:30)
accuracy <- c(5:30)

for(k in 1:length(knn.seq)) {
  fit <- CurveFinder(xy |> as.matrix(), knn=knn.seq[k])
  accuracy[k] <- cor(fit$xyt$t, t)^2
}

df <- 2:50
accuracy.pc <- df

for(k in 1:length(df)) {
  fit.ps <- principal_curve(as.matrix(xy),df=df[k])
  accuracy.pc[k] <- cor(fit.ps$lambda, t)^2
  print(accuracy.pc[k])
}


dat <- data.frame(x=knn.seq, y=accuracy)
p.knn <- ggplot(data=dat,aes(x=x,y=y)) +
  geom_point() + geom_line() +
  theme_bw() +
  ylab("r^2 with true coordinate") +
  xlab("KNN") + ggtitle("CurveFinder") +
  ylim(0,1) #+ coord_fixed()

dat <- data.frame(x=df, y=accuracy.pc)
p.df <- ggplot(data=dat,aes(x=x,y=y)) +
  geom_point() + geom_line() +
  theme_bw() +
  #ylab("r^2 with true coordinate") +
  xlab("df") + ggtitle("Principal curves") +
  ylim(0,1) #+ coord_fixed()


ggarrange(ggarrange(p, p.cp, p.p2, p.p3, nrow=1, ncol=4),
          ggarrange(p.pc.df, p.pc.df25 ,p.pc.df50 , nrow=1),
          ggarrange(p.knn, p.df,nrow=1),
          nrow=3, ncol=1)

