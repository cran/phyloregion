# [April 25, 2023] Changes in `phyloregion` version 1.0.8.

* Prepared and updated package for CRAN.

* Revised all spatial functions to use functions in `terra` throughout.

* Discontinued use of `raster`, `rgeos` and `sp` packages.

* pd and beta diversity functions now accepts a community matrix of class matrix
  
  or an `phyloseq` object as input.


# [April 30, 2021] Changes in `phyloregion` version 1.0.5

* Revised `plot_swatch` to use `col` instead of palette argument.

* Added new function, `sdm`, for species distribution modelling.

* Improved functions for transforming converting input data.


# [July 19, 2020] Changes in `phyloregion` version 1.0.3

* Removed `plot_evoldistinct` because the functionality is replaced 
with `plot_swatch`.

* `plot.phyloregion` and the generic `plot` function are now exported and both 
produce the same result. 

* Changed title to "`phyloregion`: R package for biogeographic 
regionalization and macroecology"

* Improvements in the vignettes with several examples.


# [March 30, 2020] `phyloregion` version 1.0.0

* This is the first CRAN release of phyloregion.
