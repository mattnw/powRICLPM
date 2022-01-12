
<!-- README.md is generated from README.Rmd. Please edit that file -->

# powRICLPM

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/powRICLPM)](https://CRAN.R-project.org/package=powRICLPM)
<!-- badges: end -->

`powRICLPM` is an `R` package that performs a power analysis for the
random intercept cross-lagged panel model (RI-CLPM) in a simple and
use-friendly way. It has three main functionalities:

1.  Perform an [à priori power analysis]() to get a preliminary sample
    size recommendation for detecting an effect of interest with the
    desired power level.
2.  Perform a [post hoc power analysis]() to investigate the performance
    of the RI-CLPM (w.r.t. a parameter of interest) under various
    conditions (i.e., under various sample sizes and number of repeated
    measures).
3.  Create Mplus syntax for performing a power analysis for the RI-CLPM
    using Mplus.

User guides for the above functionalities can be found in vignettes
under the ‘Articles’ tab. Technical details on the implementation of the
power analysis can be found in LINK NAAR PAPER.

## Installation

You can install the development version of `powRICLPM` from GitHub with:

``` r
install.packages("devtools")
devtools::install_github("jeroendmulder/powRICLPM")
```

## Documentation

Every user-facing function in the package is documented, and the
documentation can be accessed by running `?function_name` in the R
console, e.g., `?powRICLPM`. Furthermore, there are three main vignettes
(accessable via the ‘Articles’ tab), describing the three main
functionalities of this package, as described above.

## Citing `powRICLPM`

You can cite the R-package with the following citation:

> Mulder, J.D., (n.d.). *Performing power analysis for the RI-CLPM*

## Contact

If you have ideas, comments, or issues you would like to raise, please
get in touch.

-   Issues and idea can be raised on GitHub via
    <https://github.com/jeroendmulder/powRICLPM>
-   Pull request can be raised on GitHub via
    <https://github.com/jeroendmulder/powRICLPM/pulls>
