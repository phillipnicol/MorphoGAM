makeSwissRoll <- function(n=10^3,sd=0.03) {
  t <- runif(n=n, min=0, max=1)

  r <- 0.5 + 0.5*t
  xy <- data.frame(x=r*cos(10*t)+rnorm(n=n,sd=sd),
                   y=r*sin(10*t) + rnorm(n=n,sd=sd))
  return(as.matrix(xy))
}

rcoord_scale <- function(x) {
  min_val <- min(x)
  max_val <- max(x)

  # If all values are the same
  if (min_val == max_val) {
    return(rep(0.5, length(x)))
  }

  out <- numeric(length(x))

  if (min_val < 0 & max_val > 0) {
    # Case: spans both negative and positive
    neg <- x <= 0
    pos <- x >= 0

    out[neg] <- 0.5 * (x[neg] - min_val) / (0 - min_val)
    out[pos] <- 0.5 + 0.5 * (x[pos] - 0) / (max_val - 0)

  } else if (max_val <= 0) {
    # Case: all nonpositive (0 is outside range, above max)
    out <- (x - min_val) / (max_val - min_val)  # [0,1]

  } else if (min_val >= 0) {
    # Case: all nonnegative (0 is outside range, below min)
    out <- (x - min_val) / (max_val - min_val)  # [0,1]

  }

  return(out)
}
