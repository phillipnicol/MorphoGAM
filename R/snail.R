
set.seed(1)
n <- 10^3

t <- runif(n=n, min=0, max=1)

r <- 0.5 + 0.5*t
xy <- data.frame(x=r*cos(10*t)+rnorm(n=n,sd=0.01),
                 y=r*sin(10*t) + rnorm(n=n,sd=0.01))

p <- xy |> ggplot(aes(x=x,y=y)) + geom_point(size=0.5) +
  theme_bw() +
  ggtitle("Snail distribution")

fit.ps <- principal_curve(as.matrix(xy))
fit.ps$lambda <- fit.ps$lambda/max(fit.ps$lambda)



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


p2 <- data.frame(x=xy[,1],y=xy[,2],
                 color=fit.ps$lambda) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),colors=c("blue","grey90","red")) +
  theme_bw() +
  labs(color="Pseudotime") +
  ggtitle("Principal curves pseudotime")


rownames(xy) <- 1:n
fit <- CurveSearcher(as.matrix(xy),knn=10,cutoff=10)
p3 <- fit$plot + ggtitle("ISOMAP curve")

p4 <- data.frame(x=xy[,1],y=xy[,2],
                 color=fit$t) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point() +
  scale_color_gradientn(values=c(0,0.5,1),colors=c("blue","grey90","red")) +
  theme_bw() +
  labs(color="Trajcetory") +
  ggtitle("ISOMAP-based Trajcetory")


my.fit <- scms(Y0=t(as.matrix(xy)),X=t(as.matrix(xy)), d=1,h=0.1,minCos=0.0001)
p5 <- data.frame(x=xy[,1],y=xy[,2]) |>
  ggplot(aes(x=x,y=y)) +
  geom_point() +
  geom_point(data=data.frame(x=my.fit$Y[1,], y=my.fit$Y[2,]),col="red")


