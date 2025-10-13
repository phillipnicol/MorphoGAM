## Version 1.0.0 December 19, 2025

Initial release.

## Version 1.1.0 October 13, 2025

This update adds new functionality as well as improving the existing code.

**New functionality:**

-   Added functional PCA (based on SVD) to the `MorphoGAM()` function. These can be visualized using `plotFPCloading()` .

**Changes to existing functionality:**

-   `CurveFinderInteractive()` has been completely rewritten.
-   Option to include known list of knots, via `knots.t` in `MorphoGAM()`.
