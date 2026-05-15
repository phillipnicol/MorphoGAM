CurveFinderInteractive2 <- function(
    xy,
    loop = FALSE,
    image_path = NULL,
    flip_y = TRUE,
    flip_x = FALSE
) {
  stopifnot(is.matrix(xy) || is.data.frame(xy))
  xy <- as.matrix(xy)

  if (ncol(xy) != 2) {
    stop("`xy` must have exactly two columns (x, y).")
  }

  colnames(xy) <- c("x", "y")

  ui <- shiny::fluidPage(
    shiny::tags$h4("Click points to sketch a curve; press Smooth when done"),
    shiny::actionButton(inputId = "smooth", label = "Smooth"),
    shiny::plotOutput("plot", click = "plot_click", height = "650px"),
    shiny::br(),
    shiny::actionButton("clear", "Clear clicks")
  )

  server <- function(input, output, session) {
    vals <- shiny::reactiveValues(
      clicks = matrix(numeric(0), ncol = 2)
    )

    img_obj <- NULL

    if (!is.null(image_path)) {
      img_obj <- read_histology_image(image_path)

      if (flip_x) {
        xy[, 1] <- img_obj$width - xy[, 1]

        # flip image horizontally
        img_obj$img <- img_obj$img[, ncol(img_obj$img):1, , drop = FALSE]
      }

      if (flip_y) {
        xy[, 2] <- img_obj$height - xy[, 2]

        # flip image vertically
        img_obj$img <- img_obj$img[nrow(img_obj$img):1, , , drop = FALSE]
      }
    }

    output$plot <- shiny::renderPlot({
      if (!is.null(img_obj)) {
        plot(
          NA,
          xlim = c(0, img_obj$width),
          ylim = c(0, img_obj$height),
          xlab = "x (px)",
          ylab = "y (px)",
          asp = 1
        )

        rasterImage(
          img_obj$img,
          xleft = 0,
          ybottom = 0,
          xright = img_obj$width,
          ytop = img_obj$height
        )

        points(
          xy[, 1],
          xy[, 2],
          col = "grey60",
          pch = 16,
          cex = 0.6
        )

      } else {
        plot(
          xy,
          xlab = "x",
          ylab = "y",
          col = "grey60",
          asp = 1
        )
      }

      if (nrow(vals$clicks) > 0) {
        points(vals$clicks, pch = 19, col = "red")

        if (nrow(vals$clicks) > 1) {
          lines(
            vals$clicks[, 1],
            vals$clicks[, 2],
            col = "dodgerblue3",
            lwd = 2
          )
        }
      }
    })

    shiny::observeEvent(input$plot_click, {
      coords <- c(input$plot_click$x, input$plot_click$y)

      if (all(is.finite(coords))) {
        vals$clicks <- rbind(
          vals$clicks,
          matrix(coords, nrow = 1)
        )
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$clear, {
      vals$clicks <- matrix(numeric(0), ncol = 2)
    })

    shiny::observeEvent(input$smooth, {
      if (nrow(vals$clicks) < 3) {
        shiny::showNotification(
          "Please click at least 3 points before smoothing.",
          type = "warning"
        )
        return(NULL)
      }

      print(vals$clicks)

      res <- interactiveCurve2(
        clicks = vals$clicks,
        loop = loop,
        xy = xy,
        img = img_obj
      )

      # Return results back to original, unflipped coordinate system
      if (!is.null(img_obj)) {
        if (flip_x) {
          res$xyt$x <- img_obj$width - res$xyt$x
          res$xyt$f1 <- img_obj$width - res$xyt$f1
        }

        if (flip_y) {
          res$xyt$y <- img_obj$height - res$xyt$y
          res$xyt$f2 <- img_obj$height - res$xyt$f2
        }
      }

      shiny::stopApp(res)
    })
  }

  app <- shiny::shinyApp(ui = ui, server = server)

  if ("runGadget" %in% getNamespaceExports("shiny")) {
    fit <- shiny::runGadget(app, stopOnCancel = FALSE)
  } else {
    fit <- shiny::runApp(app)
  }

  return(fit)
}



read_histology_image <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("tif", "tiff")) {
    if (!requireNamespace("tiff", quietly = TRUE)) {
      stop("Install the 'tiff' package to read .tif/.tiff files.")
    }
    img <- tiff::readTIFF(path, native = FALSE)  # array [h,w,channels] or [h,w]
  } else if (ext == "png") {
    if (!requireNamespace("png", quietly = TRUE)) {
      stop("Install the 'png' package to read .png files.")
    }
    img <- png::readPNG(path)
  } else if (ext %in% c("jpg", "jpeg")) {
    if (!requireNamespace("jpeg", quietly = TRUE)) {
      stop("Install the 'jpeg' package to read .jpg/.jpeg files.")
    }
    img <- jpeg::readJPEG(path)
  } else {
    stop("Unsupported image type: ", ext)
  }

  # ensure 3 channels (RGB) if grayscale
  if (length(dim(img)) == 2) {
    img <- array(rep(img, 3), dim = c(dim(img), 3))
  }

  list(img = img, height = dim(img)[1], width = dim(img)[2])
}










interactiveCurve2 <- function(clicks, loop, xy, img=img) {
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
  df.line <- data.frame(x=ft[,1], y=ft[,2], color=my.t)
  h <- img$height
  w <- img$width
  p <- ggplot() +
    annotation_raster(
      raster = as.raster(img$img),
      xmin = 0, xmax = w,
      ymin = 0, ymax = h
    ) +
    geom_point(data = df.new, aes(x = x, y = y),
               col = "grey", alpha = 0.5, size = 0.6) +
    geom_path(data = df.line, aes(x = x, y = y, color = color),
              linewidth = 1) +
    scale_color_gradient(low = "navyblue", high = "firebrick1") +
    labs(color = "t") +
    coord_fixed(xlim = c(0, w), ylim = c(0, h), expand = FALSE) +
    theme_void()  # looks nicer for images; or theme_bw()

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



