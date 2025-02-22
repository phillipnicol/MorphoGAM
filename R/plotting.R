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
  Y.sub <- Y[genes,] |> as.matrix()
  offset <- Matrix::colSums(Y)
  nmed <- median(offset)
  Y.sub <- sweep(Y.sub, MARGIN = 2, FUN = "/", STATS = Matrix::colSums(Y)/nmed)
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
  return(p)
}




#' @export
#'
#' @title plotFPCloading
#'
plotFPCloading <- function(mgam_object,
                           curve.fit,
                           L=1,
                           num_genes=5,
                           type="t") {

  if(type == "t") {
    top5 <- order(abs(mgam$fpca.t$u[,L]), decreasing = TRUE)[1:num_genes]

    if(sum(mgam$fxs.t[top5,]%*%mgam$fpca.t$v[,L]) < 0) {
      mgam$fpca.t$v[,L] <- -1*mgam$fpca.t$v[,L]
    }

    mat <- matrix(0, nrow=ncol(mgam$fxs.t), ncol=num_genes + 1)
    mat[,1] <- mgam$fpca.t$d[L]*mgam$fpca.t$v[,L]
    colnames(mat) <- rep("Eigenfn", num_genes+1)
    for(i in 1:num_genes) {
      mat[,i+1] <- mgam$fxs.t[top5[i],]
      colnames(mat)[i+1] <- rownames(mgam$fxs.t)[top5[i]]
    }

    mat <- mat[,-1]
    df <- as.data.frame(mat)
    df$t <- curve.fit$xyt$t
    df <- df |> pivot_longer(cols=-t)

    p <- ggplot(data=df,aes(x=t,y=value, group=name)) +
      geom_line()

    # Generate a bright palette that avoids yellow
    #num_genes <- length(unique(df$name))
    base_colors <- brewer.pal(9, "Set1")[-6]  # Exclude the 6th color (yellow)

    # Extend the palette dynamically if more colors are needed
    if (num_genes <= length(base_colors)) {
      custom_colors <- base_colors[1:num_genes]
    } else {
      custom_colors <- colorRampPalette(base_colors)(num_genes)
    }

    # Ensure "Eigenfn" is black
    #names(custom_colors) <- colnames(mat)[-1]
    #custom_colors["Eigenfn"] <- "black"
    names(custom_colors) <- colnames(mat)

    label_data <- df %>%
      group_by(name) %>%
      slice_max(value, with_ties=FALSE,n=1)

    #print(label_data)

    # Plot
    p <- ggplot(data = df, aes(x = t, y = value,
                               group = name,
                               color = name)) +
      geom_line(aes(size = ifelse(name == "Eigenfn", 2, 0.75)),
                show.legend = FALSE) +  # Conditional size
      scale_color_manual(values = custom_colors) +  # Custom color palette
      scale_size_identity() +  # Use size as is, without scaling
      theme_bw() +
      geom_text_repel(data = label_data,
                      aes(label = name),
                      size = 4,
                      nudge_x = 0.1,
                      hjust = 0,
                      show.legend = FALSE,
                      max.overlaps = Inf) +
      ylab("Log FC from baseline") +
      geom_abline(slope=0,intercept=0,color="grey",linetype="dashed")
    xlab("t")

    return(p)
  } else{
    top5 <- order(abs(mgam$fpca.r$u[,L]), decreasing = TRUE)[1:num_genes]

    if(sum(mgam$fxs.r[top5,]%*%mgam$fpca.r$v[,L]) < 0) {
      mgam$fpca.r$v[,L] <- -1*mgam$fpca.r$v[,L]
    }

    mat <- matrix(0, nrow=ncol(mgam$fxs.r), ncol=num_genes + 1)
    mat[,1] <- mgam$fpca.r$d[L]*mgam$fpca.r$v[,L]
    colnames(mat) <- rep("Eigenfn", num_genes+1)
    for(i in 1:num_genes) {
      mat[,i+1] <- mgam$fxs.r[top5[i],]
      colnames(mat)[i+1] <- rownames(mgam$fxs.r)[top5[i]]
    }

    mat <- mat[,-1]
    df <- as.data.frame(mat)
    df$r <- curve.fit$xyt$r
    df <- df |> pivot_longer(cols=-r)

    p <- ggplot(data=df,aes(x=r,y=value, group=name)) +
      geom_line()

    # Generate a bright palette that avoids yellow
    #num_genes <- length(unique(df$name))
    base_colors <- brewer.pal(9, "Set1")[-6]  # Exclude the 6th color (yellow)

    # Extend the palette dynamically if more colors are needed
    if (num_genes <= length(base_colors)) {
      custom_colors <- base_colors[1:num_genes]
    } else {
      custom_colors <- colorRampPalette(base_colors)(num_genes)
    }

    # Ensure "Eigenfn" is black
    #names(custom_colors) <- colnames(mat)[-1]
    #custom_colors["Eigenfn"] <- "black"
    names(custom_colors) <- colnames(mat)

    label_data <- df %>%
      group_by(name) %>%
      slice_max(value, with_ties=FALSE,n=1)


    # Plot
    p <- ggplot(data = df, aes(x = r, y = value,
                               group = name,
                               color = name)) +
      geom_line(aes(size = ifelse(name == "Eigenfn", 2, 0.75)),
                show.legend = FALSE) +  # Conditional size
      scale_color_manual(values = custom_colors) +  # Custom color palette
      scale_size_identity() +  # Use size as is, without scaling
      theme_bw() +
      geom_text_repel(data = label_data,
                      aes(label = name),
                      size = 4,
                      nudge_x = 0.1,
                      hjust = 0,
                      show.legend = FALSE,
                      max.overlaps = Inf) +
      ylab("Log FC from baseline") +
      xlab("r")

    return(p)
  }
}


