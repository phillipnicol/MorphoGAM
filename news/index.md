# Changelog

## Version 1.2.0 July 11, 2026

This update adds new functionality as well as improving the existing
code.

**New functionality:**

- Add automatic selection of `k` in
  [`CurveFinder()`](https://phillipnicol.github.io/MorphoGAM/reference/CurveFinder.md)
  using `knn = "auto"`.

- A new function
  [`CurveFinderInteractive_img()`](https://phillipnicol.github.io/MorphoGAM/reference/CurveFinderInteractive_img.md)
  that allows a curve to be drawn on an image (such as a histology
  image).

**Changes to existing functionality:**

- In
  [`MorphoGAM()`](https://phillipnicol.github.io/MorphoGAM/reference/MorphoGAM.md)
  the default basis is now cubic regression splines `bs = "cr"`.

## Version 1.1.1 November 20, 2025

This update corrects some bugs in
[`CurveFinderInteractive()`](https://phillipnicol.github.io/MorphoGAM/reference/CurveFinderInteractive.md).

## Version 1.1.0 October 13, 2025

This update adds new functionality as well as improving the existing
code.

**New functionality:**

- Added functional PCA (based on SVD) to the
  [`MorphoGAM()`](https://phillipnicol.github.io/MorphoGAM/reference/MorphoGAM.md)
  function. These can be visualized using
  [`plotFPCloading()`](https://phillipnicol.github.io/MorphoGAM/reference/plotFPCloading.md)
  .
- Added code in
  [`CurveFinder()`](https://phillipnicol.github.io/MorphoGAM/reference/CurveFinder.md)
  to detect whether tissue structure is open or closed curve.

**Changes to existing functionality:**

- [`CurveFinderInteractive()`](https://phillipnicol.github.io/MorphoGAM/reference/CurveFinderInteractive.md)
  has been completely rewritten.
- Option to include known list of knots, via `knots.t` in
  [`MorphoGAM()`](https://phillipnicol.github.io/MorphoGAM/reference/MorphoGAM.md).

## Version 1.0.0 December 19, 2024

Initial release.
