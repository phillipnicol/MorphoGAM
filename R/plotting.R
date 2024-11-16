#' @export
#'
#' @title plotGAMestimates
#'
plotGAMestimates <- function(Y,
                             genes,
                             mgam_object,
                             curve_fit,
                             type="t",
                             nrow=1) {
  Y.sub <- Y[genes,]
  offset <- colSums(Y)
  nmed <- median(offset)
  Y.sub <- sweep(Y.sub, MARGIN = 2, FUN = "/", STATS = colSums(Y)/nmed)
  if(type == "t") {
    colnames(Y.sub) <- curve_fit$xyt$t
    beta_g0 <- mgam_object$results[genes,]$intercept
    fxs.sub <- nmed*exp(beta_g0 + mgam_object$fxs.t[genes,])
    colnames(fxs.sub) <- colnames(Y.sub)

    df.gg <- reshape2::melt(t(Y.sub))
    df.gg$fx <- reshape2::melt(t(fxs.sub))$value

    p <- ggplot(data=df.gg,aes(x=Var1, y=value)) +
      geom_point(size=0.25) +
      geom_line(aes(x=Var1, y=fx), color="red") +
      facet_wrap(~Var2,scales="free_y",nrow=nrow) +
      theme_bw() +
      ylab("Depth-normalized count") +
      xlab("t")
  } else if(type == "r") {
    colnames(Y.sub) <- curve_fit$xyt$r
    beta_g0 <- mgam_object$results[genes,]$intercept
    fxs.sub <- nmed*exp(beta_g0 + mgam_object$fxs.r[genes,])
    colnames(fxs.sub) <- colnames(Y.sub)

    df.gg <- reshape2::melt(t(Y.sub))
    df.gg$fx <- reshape2::melt(t(fxs.sub))$value
    print(head(df.gg))
    p <- ggplot(data=df.gg,aes(x=Var1, y=value)) +
      geom_point(size=0.25) +
      geom_line(aes(x=Var1, y=fx), color="red") +
      facet_wrap(~Var2,scales="free_y",nrow=nrow) +
      theme_bw() +
      ylab("Depth-normalized count") +
      xlab("r")
  } else{
    stop("Type must be t or r")
  }

  print(p)
  return(p)
}
