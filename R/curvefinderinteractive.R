#' @import shiny
#' @import princurve
#' @import ggplot2
#'
#' @title Interactively sketch and smooth a curve through 2D points
#'
#' @description
#' Launches a small Shiny app where you can click points to sketch a polyline
#' through 2D data \code{xy}. When you press **Smooth**, the clicked polyline is
#' smoothed with \pkg{mgcv} and then used to project the data onto the curve
#' (via \pkg{princurve}), returning coordinates \eqn{t \in [0,1]} along the
#' curve and a residual-like coordinate \eqn{r} orthogonal to it.
#'
#' @details
#' - Input \code{xy} must be 2 columns corresponding to \code{x} and \code{y}.
#' - Click at least 3 points to define a path; use **Clear clicks** to restart.
#' - Smoothing uses \code{mgcv::gam} with a cyclic cubic spline (\code{bs = "cc"})
#'   if \code{loop = TRUE}, otherwise a cubic regression spline (\code{bs = "cr"}).
#' - After smoothing, points are projected to the smoothed curve using
#'   \code{princurve::project_to_curve}, and a set of outputs (including plots)
#'   is returned and the app closes.
#'
#' @param xy A numeric matrix or data.frame with exactly two columns
#'   (interpreted as \code{x}, \code{y}). Row order is treated as the set of
#'   points to visualize/project.
#' @param loop Logical; if \code{TRUE}, the smoothed curve uses a cyclic basis
#'   (appropriate for closed loops). If \code{FALSE}, uses a non-cyclic basis.
#'
#' @return A list with the elements produced by \code{interactiveCurve()}:
#' \itemize{
#'   \item \code{xyt}: a data.frame with columns \code{x}, \code{y},
#'         \code{t} (curve coordinate scaled to \code{[0,1]}),
#'         \code{r} (orthogonal residual-like coordinate),
#'         \code{f1}, \code{f2} (fitted \code{x(t)} and \code{y(t)} values).
#'   \item \code{curve.plot}: a \pkg{ggplot2} object showing the data and the
#'         smoothed curve colored by \code{t}.
#'   \item \code{coordinate.plot}: a \pkg{ggplot2} scatter plot colored by \code{t}.
#'   \item \code{residuals.plot}: a \pkg{ggplot2} scatter plot colored by \code{r}.
#' }
#'
CurveFinderInteractive <- function(xy, loop = FALSE) {
  stopifnot(is.matrix(xy) || is.data.frame(xy))
  xy <- as.matrix(xy)
  if (ncol(xy) != 2) stop("`xy` must have exactly two columns (x, y).")
  colnames(xy) <- c("x", "y")

  ui <- shiny::fluidPage(
    shiny::tags$h4("Click points to sketch a curve; press Smooth when done"),
    shiny::actionButton(inputId = "smooth", label = "Smooth"),
    shiny::plotOutput("plot", click = "plot_click", height = "500px"),
    shiny::br(),
    shiny::actionButton("clear", "Clear clicks")
  )

  server <- function(input, output, session) {
    vals <- shiny::reactiveValues(clicks = matrix(numeric(0), ncol = 2))

    output$plot <- shiny::renderPlot({
      plot(xy, xlab = "x", ylab = "y", col = "grey60")
      if (nrow(vals$clicks) > 0) {
        points(vals$clicks, pch = 19, col = "red")
        if (nrow(vals$clicks) > 1) {
          lines(vals$clicks[,1], vals$clicks[,2], col = "blue", lwd = 2)
        }
      }
    })

    shiny::observeEvent(input$plot_click, {
      coords <- c(input$plot_click$x, input$plot_click$y)
      if (all(is.finite(coords))) {
        vals$clicks <- rbind(vals$clicks, matrix(coords, nrow = 1))
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$clear, {
      vals$clicks <- matrix(numeric(0), ncol = 2)
    })

    shiny::observeEvent(input$smooth, {
      if (nrow(vals$clicks) < 3) {
        shiny::showNotification("Please click at least 3 points before smoothing.", type = "warning")
        return(NULL)
      }
      res <- interactiveCurve(clicks = vals$clicks, loop = loop, data = xy)
      shiny::stopApp(res)
    })
  }

  app <- shiny::shinyApp(ui = ui, server = server)
  # If available, run as gadget; otherwise run as normal app.
  if ("runGadget" %in% getNamespaceExports("shiny")) {
    fit <- shiny::runGadget(app, stopOnCancel = FALSE)
  } else {
    fit <- shiny::runApp(app)
  }
  return(fit)
}


interactiveCurve <- function(clicks, loop, data) {
  print("Running smoother. Do not press stop.")
  #my.clicks <- interactiveCurve(xy)
  #clicks <- interactiveCurve(xy)
  basis <- ifelse(loop, "cc", "cr")
  t <- seq(0, 1, length.out=nrow(clicks))
  fitx <- mgcv::gam(clicks[,1]~s(t,bs=basis, k=nrow(clicks)))
  fity <- mgcv::gam(clicks[,2]~s(t,bs=basis, k=nrow(clicks)))
  my.t <- seq(0,1,by=10^{-4})
  predx <- predict(fitx,newdata=list(t=my.t))
  predy <- predict(fity,newdata=list(t=my.t))
  ft <- data.frame(x=predx,y=predy)

  proj <- princurve::project_to_curve(xy,s=as.matrix(ft))
  t <- as.numeric(proj$lambda/max(proj$lambda))

  # cumulative arc length along the dense curve 'ft'
  arc <- c(0, cumsum(sqrt(diff(ft$x)^2 + diff(ft$y)^2)))
  arc <- arc / max(arc)                     # normalize to [0,1] like proj$lambda/max

  # interpolate lambda -> t on the dense grid
  t_at_proj <- approx(x = arc, y = my.t,
                      xout = proj$lambda / max(proj$lambda),
                      rule = 2)$y

  #fitx <- mgcv::gam(xy[,1]~s(t,bs=basis, k=nrow(clicks)))
  #fity <- mgcv::gam(xy[,2]~s(t,bs=basis, k=nrow(clicks)))

  r <- orthogonal_path_interactive(fitx,fity,t_at_proj,xy)

  df.new <- data.frame(x=xy[,1],
                       y=xy[,2])

  p <- ggplot(data=df.new,aes(x=x,y=y)) + geom_point(col="grey",
                                                     alpha=0.5)

  df.line <- data.frame(x=ft[,1], y=ft[,2], color=my.t)
  p <- p + geom_path(data=df.line,aes(x=x,y=y,color=color),
                     linewidth=1)
  p <- p + scale_color_gradient(low="navyblue",
                                high="firebrick1")
  p <- p + theme_bw() + labs(color="t")

  xyt <- data.frame(x=xy[,1],y=xy[,2],t=t, r=r,
                    f1 = predict(fitx, newdata=list(t=t)),
                    f2=predict(fitx, newdata=list(t=t)))


  p2 <- data.frame(x=xy[,1],y=xy[,2],color=t) |>
    ggplot(aes(x=x,y=y,color=color)) + geom_point() +
    scale_color_gradientn(values=c(0,0.5,1),
                          colors=c("blue","grey90", "red"))+
    theme_bw() +
    ggtitle("Coordinate") + labs(color="t")

  p3 <- data.frame(x=xy[,1],y=xy[,2],color=r) |>
    ggplot(aes(x=x,y=y,color=color)) + geom_point() +
    scale_color_gradientn(values=c(0,0.5,1),
                          colors=c("blue","grey90", "red"))+
    theme_bw() +
    ggtitle("Curve residuals") + labs(color="r")

  out <- list()
  out$xyt <- xyt
  out$curve.plot <- p
  out$coordinate.plot <- p2
  out$residuals.plot <- p3
  print("Done!")
  return(out)
}



orthogonal_path_interactive <- function(fitx,fity,t,xy) {
  f2x <- gratia::derivatives(fitx,order=1,data=data.frame(t=t))
  f2y <- gratia::derivatives(fity,order=1,data=data.frame(t=t))
  t2 <- t

  x.fit <- predict(fitx, newdata=data.frame(t=t))
  y.fit <- predict(fity, newdata=data.frame(t=t))

  residual.x <- xy[,1] - x.fit
  residual.y <- xy[,2] - y.fit

  for(i in 1:length(t)) {
    e <- c(residual.x[i], residual.y[i])
    Rf1 <- c(-f2y$.derivative[i], f2x$.derivative[i])
    sign <- ifelse(sum(e*Rf1) > 0, 1, -1)
    t2[i] <- sign*sqrt(sum(e^2))
    #t2[i] <- sqrt(sum(e^2))
  }

  #t2 <- (t2 - min(t2))/(max(t2) - min(t2))
  t2 <- rcoord_scale(t2)
  return(t2)
}


