
get.clicks <- function(xy) {
  interactiveCurve(xy)
}

CurveFinderInteractive <- function(clicks, loop=FALSE) {
  #my.clicks <- interactiveCurve(xy)

  basis <- ifelse(loop, "cc", "cr")
  print("BACK")
  t <- seq(0, 1, length.out=nrow(clicks))
  print("BACK")
  fitx <- mgcv::gam(clicks[,1]~s(t,bs=basis, k=nrow(clicks)))
  fity <- mgcv::gam(clicks[,2]~s(t,bs=basis, k=nrow(clicks)))
  print("BACK")
  my.t <- seq(0,1,by=10^{-4})
  print("BACK")
  predx <- predict(fitx,newdata=list(t=my.t))
  predy <- predict(fity,newdata=list(t=my.t))
  print("BACK")
  ft <- data.frame(x=predx,y=predy)

  print("BACK")
  proj <- project_to_curve(xy,s=as.matrix(ft))
  t <- as.numeric(proj$lambda/max(proj$lambda))
  print(t)
  print(xy)

  fitx <- mgcv::gam(xy[,1]~s(t,bs=basis, k=nrow(clicks)))
  fity <- mgcv::gam(xy[,2]~s(t,bs=basis, k=nrow(clicks)))

  print("BACK")
  r <- orthogonal_path(fitx,fity,t)

  df.new <- data.frame(x=xy[,1],
                       y=xy[,2])

  p <- ggplot(data=df.new,aes(x=x,y=y)) + geom_point(col="grey",
                                                     alpha=0.5)

  df.line <- data.frame(x=ft[,1], y=ft[,2], color=my.t)
  p <- p + geom_path(data=df.line,aes(x=x,y=y,color=color),
                     linewidth=1)
  p <- p + scale_color_gradient(low="navyblue",
                                high="firebrick1")
  p <- p + theme_bw()

  print("BACK")
  xyt <- data.frame(x=xy[,1],y=xy[,2],t=t, r=r,
                    f1 = fitted(fitx),
                    f2=fitted(fity))


  p2 <- data.frame(x=xy[,1],y=xy[,2],color=t) |>
    ggplot(aes(x=x,y=y,color=color)) + geom_point() +
    scale_color_gradientn(values=c(0,0.5,1),
                          colors=c("blue","grey90", "red"))+
    theme_bw() +
    ggtitle("Coordinate")

  p3 <- data.frame(x=xy[,1],y=xy[,2],color=r) |>
    ggplot(aes(x=x,y=y,color=color)) + geom_point() +
    scale_color_gradientn(values=c(0,0.5,1),
                          colors=c("blue","grey90", "red"))+
    theme_bw() +
    ggtitle("Curve residuals")

  out <- list()
  out$xyt <- xyt
  out$curve.plot <- p
  out$coordinate.plot <- p2
  out$residuals.plot <- p3
  return(out)
}


interactiveCurve <- function(xy) {
  ## Written with Claude
  # Define UI
  ui <- fluidPage(
    plotOutput("plot", click = "plot_click"),
    verbatimTextOutput("clicks")
  )

  my.clicks <- rep(NA, nrow=0, ncol=2)

  # Define server logic
  server <- function(input, output, session) {
    # Initialize a reactive value to store clicked points
    vals <- reactiveValues(
      clicks = NULL,
      prev_click = NULL
    )

    # Generate plot
    output$plot <- renderPlot({
      plot(xy, xlab = "x", ylab = "y")

      # Add clicked points to the plot
      if (!is.null(vals$clicks)) {
        points(vals$clicks, pch = 19, col = "red")

        # Draw line segments between consecutive clicked points
        if (nrow(vals$clicks) > 1) {
          for (i in 2:nrow(vals$clicks)) {
            x1 <- vals$clicks[i - 1, 1]
            y1 <- vals$clicks[i - 1, 2]
            x2 <- vals$clicks[i, 1]
            y2 <- vals$clicks[i, 2]
            lines(c(x1, x2), c(y1, y2), col = "blue")
          }
        }
      }
    })

    # Store clicked points
    observeEvent(input$plot_click, {
      coords <- input$plot_click
      # Update the clicked points with the new click
      vals$clicks <- rbind(vals$clicks, coords)

      # Store the previous click
      vals$prev_click <- coords
      my.clicks <<- rbind(my.clicks, coords[1:2])
    })

    session$onSessionEnded(function() {
      stopApp()
    })
  }
  print("Start")

  # Run the app
  runApp(shinyApp(ui = ui, server = server))
  print("END")
  return(matrix(unlist(my.clicks[-1,]), ncol=2))
}
