# Plot depth-normalized counts with fitted GAM overlay

For a set of `genes`, this plots depth-normalized counts across either
the coordinate `t` or the coordinate `r`, with an optional overlaid
fitted function from a mixed GAM object.

## Usage

``` r
plotGAMestimates(
  Y,
  genes,
  mgam_object,
  curve_fit,
  type = "t",
  nrow = 1,
  include.gam = TRUE
)
```

## Arguments

- Y:

  A numeric matrix-like (e.g., `matrix`, `Matrix::dgCMatrix`) of counts
  with **genes in rows** and **cells/samples in columns**. Row names
  should include the values provided in `genes`.

- genes:

  A character vector (or integer indices) of gene identifiers to plot.

- mgam_object:

  The output from running the `MorphoGAM` function.

- curve_fit:

  The output from running the `CurveFinder` (or interactive version)
  function.

- type:

  Character string; one of `"t"` or `"r"` specifying which coordinate to
  plot against.

- nrow:

  Integer; number of rows in the facet layout.

- include.gam:

  Logical; if `TRUE`, overlay the red GAM fit line.

## Value

A `ggplot2` object showing depth-normalized counts for each gene
(facets) versus `t` or `r`, with an optional fitted curve.

## Details

Columns (cells/samples) are depth-normalized by dividing by the library
size (column sum) and multiplying by the median library size across all
columns. If `type = "t"`, the x-axis uses `curve_fit$xyt$t` and fitted
values come from `mgam_object$fxs.t`. If `type = "r"`, the x-axis uses
`curve_fit$xyt$r` and fitted values come from `mgam_object$fxs.r`.
Fitted curves are constructed as \\n\_{med} \cdot \exp(\beta\_{g0} +
f_g(x))\\, where \\\beta\_{g0}\\ is the gene-specific intercept and
\\f_g(x)\\ is the smooth contribution for that gene.
