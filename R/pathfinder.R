

PathFinder <- function(xy,
                       knn=5,
                       cutoff=10) {
  xy.dist <- as.matrix(dist(xy))

  nnk <- apply(xy.dist, 1, function(x) sort(x)[knn+1])
  outlier <- which(nnk > cutoff*median(nnk))

  if(length(outlier) > 0) {
    xy.new <- xy[-outlier,]
  } else {
    xy.new <- xy
  }


}
