# Interactively sketch and smooth a curve on an optional image

Launches a Shiny app for manually sketching a curve through
two-dimensional points, optionally overlaid on a histology or spatial
image. The clicked curve is smoothed and used to assign each point a
curve coordinate \`t\` and an orthogonal residual-like coordinate \`r\`.

## Usage

``` r
CurveFinderInteractive_img(
  xy,
  loop = FALSE,
  image_path = NULL,
  flip_y = TRUE,
  flip_x = FALSE
)
```

## Arguments

- xy:

  A numeric matrix or data frame with exactly two columns containing x
  and y coordinates.

- loop:

  Logical; if \`TRUE\`, smooth the clicked curve with a cyclic basis for
  a closed loop. If \`FALSE\`, use a non-cyclic cubic regression basis.

- image_path:

  Optional path to a \`.png\`, \`.jpg\`, \`.jpeg\`, \`.tif\`, or
  \`.tiff\` image to show behind the points.

- flip_y:

  Logical; if \`TRUE\`, flip the y-axis while displaying the image and
  points. Default is \`TRUE\`.

- flip_x:

  Logical; if \`TRUE\`, flip the x-axis while displaying the image and
  points. Default is \`FALSE\`.

## Value

A list with fitted curve output:

- \`xyt\`: a data frame with \`x\`, \`y\`, \`t\`, \`r\`, \`f1\`, and
  \`f2\`.

- \`curve.plot\`: a \`ggplot2\` object showing the smoothed curve.

- \`coordinate.plot\`: a scatter plot colored by \`t\`.

- \`residuals.plot\`: a scatter plot colored by \`r\`.

## Details

If \`image_path\` is provided, PNG, JPEG, TIFF, and TIF files are
supported when the corresponding optional reader package is installed.
\`flip_x\` and \`flip_y\` transform both the displayed image and point
coordinates before curve fitting, and returned coordinates are
transformed back to the original coordinate system.
