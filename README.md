
<!-- README.md is generated from README.Rmd. Please edit that file -->
sbar - Stage-based assessments in R
===================================

<!-- badges: start -->
<!-- badges: end -->
sbar is a package that provides functions to allow the user to run various stage-based fisheries assessments in R. The stage-based assesments available are versions of the current (2020) Catch-Survey Analysis (CSA) found in the the NOAA Fish and Fisheries Toolbox and an implementation of a size-based delay-diffrenec model described in the throerteical paper by Schnute (1987). Model theory, implementation and testing are detailed in ........

Installation
------------

This package involves C++ code and the modelling framework of the package . Please ensure that you are able to compile source code and that this package is running smoothly on your device.

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

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub!
