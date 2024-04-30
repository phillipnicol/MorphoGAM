

library(STexampleData)

spe <- SlideSeqV2_mouseHPC()
spe <- spe[,!is.na(spe$celltype)]

for(ct in unique(spe$celltype)) {
  print(ct)
  ixs <- which(spe$celltype == ct)

  if(length(ixs) < 50) {
    next
  }

  xy <- spatialCoords(spe)[ixs,]

  fit <- CurveSearcher.cv(xy,nfolds=10)

  p <- fit$plot + ggtitle(paste0(ct, " ", fit$Rhat))
  ggsave(p, filename=paste0("../plots/", ct, "curve.png"))
}


ixs <- which(spe$celltype == "CA3") #subset to CA3

xy <- spatialCoords(spe)[ixs,]
Y <- counts(spe)[,ixs]
xy.dist <- as.matrix(dist(xy))

nnk <- apply(xy.dist, 1, function(x) sort(x)[3])
outlier <- which(nnk > 2*median(nnk))

if(length(outlier) > 0) {
  xy.new <- xy[-outlier,]
} else {
  xy.new <- xy
}

Y <- Y[,-outlier]


knng <- dimRed:::makeKNNgraph(x = xy.new,
                              k = knn,
                              eps = 0)


comp <- components(knng)
xy <- xy.new[comp$membership == 1,]
Y <- Y[,comp$membership == 1]

fit <- CurveSearcher(xy,knn=5,cutoff=10)

my.svg <- sparkx(count_in=Y,locus_in =xy)

p.vals <- my.svg$res_mtest |> as.data.frame() |>
  arrange(adjustedPval)

gene <- rep(0, ncol(Y))
G <- gam(gene~s(fit$t,k=10,bs="cr"), family=nb(), fit=FALSE)
X <- G$X

t <- fit$t
l.o <- log(colSums(Y))
coef <- rep(0, nrow(Y))
p.val <- rep(1,nrow(Y))
nz <- apply(Y, 1, function(x) sum(x != 0))
f <- matrix(0,nrow=nrow(Y), ncol=length(t))
for(i in 1:nrow(Y)) {
  if(nz[i] < 50) {
    next
  }
  print(i)
  my.fit <- glm(Y[i,]~t+offset(l.o), family=quasipoisson())
  p.val[i] <- summary(my.fit)$coefficients[2,4]
  coef[i] <- summary(my.fit)$coefficients[2,1]
}

df <- data.frame(p.vals=p.val, slope=coef, spark_rank = match(rownames(Y),rownames(p.vals)))
rownames(df) <- rownames(Y)

df |> filter(p.vals < 10^{-6}) |> arrange(desc(slope)) |> head(n=10)
df |> filter(p.vals < 10^{-6}) |> arrange(slope) |> head(n=10)

