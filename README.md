
<!-- README.md is generated from README.Rmd. Please edit that file -->
sbar - Stage-based assessments in R
===================================

<!-- badges: start -->
<!-- badges: end -->
sbar is a package that provides functions to allow the user to run various stage-based fisheries assessments in R. The stage-based assesments available are versions of the current (2020) Catch-Survey Analysis (CSA) found in the the NOAA Fish and Fisheries Toolbox and an implementation of a size-based delay-diffrenec model described in the throerteical paper by Schnute (1987). Model theory, implementation and testing are detailed in ........

Installation
------------

This package involves C++ code and the modelling framework of the package . Please ensure that you are able to compile source code and that this package is running smoothly on your device. The package has not been tested on the most recent version of R (>=4.0) as of yet

And sbar package can be installed from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lbatts/sbar")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(sbar)
## basic example code
```
You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

