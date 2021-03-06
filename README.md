
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
use-friendly way. Its main functionalities are:

1.  [Perform a power
    analysis](https://jeroendmulder.github.io/powRICLPM/articles/analysis.html)
    to obtain sample size recommendations (as well as other performance
    measures, such as bias, mean square error, etc.) for a desired power
    level for a specific parameter. This can be done across multiple
    experimental conditions simultaneously (i.e., across varying numbers
    of repeated measures, proportions of between-unit variance, etc.).
2.  [Summarize and visualize power analysis
    results](https://jeroendmulder.github.io/powRICLPM/reference/visualization.html).
3.  [Create Mplus
    syntax](https://jeroendmulder.github.io/powRICLPM/articles/mplus.html)
    for performing a power analysis for the RI-CLPM using Mplus.

User guides for the above functionalities can be found in vignettes
under the ‘Articles’ tab. Technical details on the implementation of the
power analysis can be found in the package’s and functions’
documentation. A rationale for the implemented power analysis strategy
here, as well as an illustrative example and extensions to extensions of
the RI-CLPM, see Mulder (forthcoming).

## Installation

You can install the development version of `powRICLPM` from GitHub with:

``` r
install.packages("devtools")
devtools::install_github("jeroendmulder/powRICLPM")
```

## Documentation

Every user-facing function in the package is documented, and the
documentation can be accessed by running `?function_name` in the R
console, e.g., `?powRICLPM`. Furthermore, there are four main vignettes
(accessible via the ‘Articles’ tab), describing functionalities and
analysis options of this package.

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
