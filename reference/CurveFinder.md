# CurveFinder: Automatically estimate curve in two-dimensional ST data

This function estimates a (possibly closed) parametric curve passing
through spatial transcriptomics data and uses this curve to define
morphologically relevant coordinates \`t\` and \`r\`.

## Usage

``` r
CurveFinder(
  xy,
  knn = 5,
  prune.outlier = NULL,
  loop = FALSE,
  knot.fx = 100,
  scale_morpho_coords = TRUE,
  max.comp.no = 5
)
```

## Arguments

- xy:

  A numeric matrix or data frame with two columns representing the x and
  y coordinates of the data points.

- knn:

  An integer specifying the number of nearest neighbors used to
  construct the KNN graph, or \`"auto"\` to choose a value from 3 to 20
  by model score. Default is 5.

- prune.outlier:

  A numeric threshold for pruning outliers based on the distance to the
  k+1 nearest neighbor. Outliers are removed if their distance exceeds
  \`prune.outlier \* median(nnk)\`. Defaults to NULL (no pruning).

- loop:

  A logical value indicating whether the curve should be treated as a
  loop (closed curve), or \`"auto"\` to detect loops from a Mapper
  graph. Default is FALSE.

- knot.fx:

  Maximum number of knots used for the GAM smooths fitted to x(t) and
  y(t). The effective value is capped at 10 percent of the number of
  points.

- scale_morpho_coords:

  A logical value indicating whether to scale the morphologically
  relevant coordinates \`t\` and \`r\` to the range \[0, 1\]. Default is
  TRUE.

- max.comp.no:

  An integer specifying the maximum number of disconnected components
  allowed in the KNN graph. Default is 5.

## Value

A list containing:

- xyt:

  A data frame with x and y coordinates, fitted curve parameters (\`t\`
  and \`r\`), and fitted values for x and y (\`f1\` and \`f2\`).

- curve.plot:

  A ggplot object showing the original data with the fitted curve
  overlaid.

- coordinate.plot:

  A ggplot object displaying the data points color-coded by their fitted
  \`t\` values.

- residuals.plot:

  A ggplot object displaying the data points color-coded by their fitted
  \`r\` values.

- model.score:

  A numeric model score used internally for \`knn = "auto"\`.

- arclength:

  The approximate arclength of the fitted curve.

- span.r:

  The estimated span of the residual coordinate.

- t_v_r_span:

  The ratio of fitted curve arclength to residual-coordinate span.

## Details

The function builds a k-nearest-neighbor graph on \`xy\`, derives an
initial one-dimensional ordering from graph distances, fits smooth
functions for the two spatial coordinates, and then computes signed
orthogonal residuals from the fitted curve.

## References

Kraemer G, Reichstein M, Mahecha MD (2018). “dimRed and
coRanking—Unifying Dimensionality Reduction in R.” \_The R Journal\_,
\*10\*(1), 342-358.
