# Plot top gene loadings from SVD Analysis

Visualize the \\L\\-th FPCA eigenfunction together with the top
`num_genes` genes most strongly associated with that component, over
either the `t` or `r` coordinate.

## Usage

``` r
plotFPCloading(mgam_object, curve.fit, L = 1, num_genes = 5, type = "t")
```

## Arguments

- mgam_object:

  The output of `MorphoGAM`

- curve.fit:

  The output of `CurveFinder`

- L:

  Integer; index of the FPCA component.

- num_genes:

  Integer; number of top-loading genes to display.

- type:

  Character string; one of `"t"` or `"r"` specifying which coordinate to
  plot against.

## Value

A `ggplot2` object showing lines for the eigenfunction and the selected
gene smooths. A bright qualitative palette (excluding yellow) is used
for gene curves; labels are placed near each curve's maximum.

## Details

The function identifies the `num_genes` genes with the largest absolute
scores in the FPCA decomposition for the requested coordinate. It may
flip sign to improve visualization.
