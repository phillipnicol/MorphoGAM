# Interactively sketch and smooth a curve through 2D points

Launches a small Shiny app where you can click points to sketch a
polyline through 2D data `xy`. When you press \*\*Smooth\*\*, the
clicked polyline is smoothed with mgcv and then used to project the data
onto the curve (via princurve), returning coordinates \\t \in \[0,1\]\\
along the curve and a residual-like coordinate \\r\\ orthogonal to it.

## Usage

``` r
CurveFinderInteractive(xy, loop = FALSE)
```

## Arguments

- xy:

  A numeric matrix or data.frame with exactly two columns (interpreted
  as `x`, `y`). Row order is treated as the set of points to
  visualize/project.

- loop:

  Logical; if `TRUE`, the smoothed curve uses a cyclic basis
  (appropriate for closed loops). If `FALSE`, uses a non-cyclic basis.

## Value

A list with the elements produced by `interactiveCurve()`:

- `xyt`: a data.frame with columns `x`, `y`, `t` (curve coordinate
  scaled to `[0,1]`), `r` (orthogonal residual-like coordinate), `f1`,
  `f2` (fitted `x(t)` and `y(t)` values).

- `curve.plot`: a ggplot2 object showing the data and the smoothed curve
  colored by `t`.

- `coordinate.plot`: a ggplot2 scatter plot colored by `t`.

- `residuals.plot`: a ggplot2 scatter plot colored by `r`.

## Details

\- Input `xy` must be 2 columns corresponding to `x` and `y`. - Click at
least 3 points to define a path; use \*\*Clear clicks\*\* to restart. -
Smoothing uses [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) with
a cyclic cubic spline (`bs = "cc"`) if `loop = TRUE`, otherwise a cubic
regression spline (`bs = "cr"`). - After smoothing, points are projected
to the smoothed curve using
[`princurve::project_to_curve`](https://rdrr.io/pkg/princurve/man/project_to_curve.html),
and a set of outputs (including plots) is returned and the app closes.
