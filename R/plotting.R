#' @export
#'
#' @title Plot depth-normalized counts with fitted GAM overlay
#'
#' @description
#' For a set of \code{genes}, this plots depth-normalized counts across either
#' the coordinate \code{t} or the coordinate \code{r}, with an
#' optional overlaid fitted function from a mixed GAM object.
#'
#' @details
#' Columns (cells/samples) are depth-normalized by dividing by the library size
#' (column sum) and multiplying by the median library size across all columns.
#' If \code{type = "t"}, the x-axis uses \code{curve_fit$xyt$t} and fitted
#' values come from \code{mgam_object$fxs.t}. If \code{type = "r"}, the x-axis
#' uses \code{curve_fit$xyt$r} and fitted values come from \code{mgam_object$fxs.r}.
#' Fitted curves are constructed as \eqn{n_{med} \cdot \exp(\beta_{g0} + f_g(x))},
#' where \eqn{\beta_{g0}} is the gene-specific intercept and \eqn{f_g(x)} is the
#' smooth contribution for that gene.
#'
#' @param Y A numeric matrix-like (e.g., \code{matrix}, \code{Matrix::dgCMatrix})
#'   of counts with \strong{genes in rows} and \strong{cells/samples in columns}.
#'   Row names should include the values provided in \code{genes}.
#' @param genes A character vector (or integer indices) of gene identifiers to plot.
#' @param mgam_object The output from running the \code{MorphoGAM} function.
#' @param curve_fit The output from running the \code{CurveFinder} (or interactive version) function.
#' @param type Character string; one of \code{"t"} or \code{"r"} specifying which
#'   coordinate to plot against.
#' @param nrow Integer; number of rows in the facet layout.
#' @param include.gam Logical; if \code{TRUE}, overlay the red GAM fit line.
#'
#' @return A \code{ggplot2} object showing depth-normalized counts for each gene
#'   (facets) versus \code{t} or \code{r}, with an optional fitted curve.
#'
plotGAMestimates <- function(Y,
                             genes,
                             mgam_object,
                             curve_fit,
                             type="t",
                             nrow=1,
                             include.gam=TRUE) {
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
      geom_point(size=0.25)

    if(include.gam) {
      p <- p + geom_line(aes(x=Var1, y=fx), color="red")
    }

    p <- p + facet_wrap(~Var2,scales="free_y",nrow=nrow) +
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
    p <- ggplot(data=df.gg,aes(x=Var1, y=value)) +
      geom_point(size=0.25)

      if(include.gam) {
        p <- p + geom_line(aes(x=Var1, y=fx), color="red")
      }

    p <- p + facet_wrap(~Var2,scales="free_y",nrow=nrow) +
      theme_bw() +
      ylab("Depth-normalized count") +
      xlab("r")
  } else{
    stop("Type must be t or r")
  }
  return(p)
}



#' @import ggplot2
#' @import RColorBrewer
#' @import ggrepel
#' @importFrom stats reorder
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr group_by slice_max
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @title Plot top gene loadings from SVD Analysis
#'
#' @description
#' Visualize the \eqn{L}-th FPCA eigenfunction together with the top
#' \code{num_genes} genes most strongly associated with that component,
#' over either the \code{t} or \code{r} coordinate.
#'
#' @details
#' The function identifies the \code{num_genes} genes with the largest absolute
#' scores in \code{mgam_object$fpca.t}. It may flip sign to improve visualization.
#'
#' @param mgam_object The output of \code{MorphoGAM}
#' @param curve.fit The output of \code{CurveFinder}
#' @param L Integer; index of the FPCA component.
#' @param num_genes Integer; number of top-loading genes to display.
#' @param type Character string; one of \code{"t"} or \code{"r"} specifying
#'   which coordinate to plot against.
#'
#' @return A \code{ggplot2} object showing lines for the eigenfunction and the
#'   selected gene smooths. A bright qualitative palette (excluding yellow) is
#'   used for gene curves; labels are placed near each curve's maximum.
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


