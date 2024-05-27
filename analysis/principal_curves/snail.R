
set.seed(1)
n <- 10^3

t <- runif(n=n, min=0, max=1)
t.true <- t

r <- 0.5 + 0.5*t
xy <- data.frame(x=r*cos(10*t)+rnorm(n=n,sd=0.01),
                 y=r*sin(10*t) + rnorm(n=n,sd=0.01))
xy <- as.matrix(xy)


fit <- CurveFinder(xy,knn=8)

p <- xy |> ggplot(aes(x=x,y=y)) + geom_point(size=0.5) +
  theme_bw() +
  ggtitle("Swiss roll")
p1 <- fit$curve.plot + guides(color="none") + ggtitle("Path (knn=8)")
p2 <- fit$coordinate.plot +guides(color="none")

library(princurve)

fit.ps <- principal_curve(as.matrix(xy))

fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)
fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p.ps.1 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("princurve (Default)") +
  guides(color="none")

fit.ps <- principal_curve(as.matrix(xy),df=50)

fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)
fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p.ps.2 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("princurve (df=50)") +
  guides(color="none")

fit.ps <- principal_curve(as.matrix(xy),df=100)

fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)
fit.ps$s <- fit.ps$s[order(fit.ps$lambda),]
df.line <- data.frame(x=fit.ps$s[,1], y=fit.ps$s[,2],
                      color=sort(fit.ps$lambda))
p.ps.3 <- data.frame(x=xy[,1], y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) + geom_point(col="grey", alpha=0.5) +
  geom_path(data=df.line,aes(x=x,y=y,color=color),
            linewidth=1) + scale_color_gradient(low="navyblue",
                                                high="firebrick1") +
  theme_bw() +
  ggtitle("princurve (df=100)") +
  guides(color="none")

knn <- seq(6,35)
res <- data.frame(x=knn, y=0)
for(k in 1:nrow(res)) {
  print(res[k,1])
  fitknn <- CurveFinder(xy,knn=res[k,1])
  res[k,2] <- cor(fitknn$xyt$t, t.true,method="spearman")^2
  print(res[k,2])
}

df <- seq(10, 100)
res2 <- data.frame(x=df, y=0)
for(i in 1:nrow(res2)) {
  fit.ps <- principal_curve(as.matrix(xy), df=res2[i,1])
  fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)
  res2[i,2] <- cor(fit.ps$lambda, t.true,method="spearman")^2
  print(res2[i,2])
}



p.3.1 <- data.frame(x=res[,1], y=res[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point() + geom_line() +
  theme_bw() + xlab("knn") + ylab("r^2 (spearman) with ground truth") +
  ylim(0,1)

p.3.2 <- data.frame(x=res2[,1], y=res2[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point() + geom_line() +
  theme_bw() + xlab("princurve df") +
  ylab("") +
  ylim(0,1)

p.big <- ggarrange(ggarrange(p,p1,p2,nrow=1,ncol=3),
                   ggarrange(p.ps.1,p.ps.2,p.ps.3, nrow=1,ncol=3),
                   ggarrange(p.3.1,p.3.2,nrow=1,ncol=2),
                   nrow=3,ncol=1,heights=c(2,2,1.5))

ggsave(p.big, file="../plots/snail.png", width=8, height=8, units="in")
