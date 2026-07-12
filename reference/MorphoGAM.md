# MorphoGAM: Apply a GAM to morphologically relevant coordinates

Fit gene-wise negative-binomial GAMs to identify expression variation
with respect to morphologically relevant coordinates estimated by
\`CurveFinder()\` or \`CurveFinderInteractive()\`.

DESCRIPTION

## Usage

``` r
MorphoGAM(
  Y,
  curve.fit,
  design,
  shrinkage = FALSE,
  min.count.per.gene = 10,
  return.fx = TRUE,
  offset = NULL,
  knots.t = NULL,
  knots.r = NULL
)
```

## Arguments

- Y:

  A numeric count matrix with genes in rows and samples or cells in
  columns.

- curve.fit:

  An object produced by \`CurveFinder()\` or
  \`CurveFinderInteractive()\`, containing fitted \`t\` and \`r\`
  coordinates for each sample.

- design:

  A formula specifying the GAM design, typically including smooth terms
  for \`t\` and \`r\` such as \`y ~ s(t) + s(r)\`. Smooth terms that do
  not specify \`bs\` use \`bs = "cr"\` by default.

- shrinkage:

  Logical; if \`TRUE\`, shrink smooth term coefficients using \`ashr\`.
  Default is \`FALSE\`.

- min.count.per.gene:

  Minimum total count required for a gene to be fit. Genes below this
  threshold are returned with null/default statistics.

- return.fx:

  Logical; if \`TRUE\`, return fitted smooth terms and FPCA summaries
  for the \`t\` and \`r\` effects.

- offset:

  Optional numeric vector of model offsets. If \`NULL\`, offsets are
  computed as the logarithm of column sums of \`Y\`.

- knots.t:

  Optional numeric vector specifying knot locations for the smooth term
  in \`t\`.

- knots.r:

  Optional numeric vector specifying knot locations for the smooth term
  in \`r\`.

## Value

A list containing:

- results:

  A data frame with one row per gene and columns \`peak.t\`,
  \`range.t\`, \`pv.t\`, \`peak.r\`, \`range.r\`, \`pv.r\`, and
  \`intercept\`.

- design:

  The model formula used after adding default cubic regression spline
  bases to unspecified smooth terms.

If \`return.fx = TRUE\`, the list also includes:

- fxs.t:

  A matrix of fitted smooth terms for \`t\` (genes by samples).

- fxs.r:

  A matrix of fitted smooth terms for \`r\` (genes by samples).

- fpca.t:

  An \`irlba\` decomposition of \`fxs.t\`, when available.

- fpca.r:

  An \`irlba\` decomposition of \`fxs.r\`, when available.

## See also

Useful links:

- <https://phillipnicol.github.io/MorphoGAM/>

## Author

Phillip B. Nicol \<philnicol740@gmail.com\>
