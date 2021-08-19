
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sbar - Stage-based assessments in R

<!-- badges: start -->
<!-- badges: end -->

sbar is a package that provides functions to allow the user to run
various stage-based fisheries assessments in R. The stage-based
assesments available are versions of the current (2020) Catch-Survey
Analysis (CSA) found in the the NOAA Fish and Fisheries Toolbox and an
implementation of a size-based delay-difference model described in the
throerteical paper by Schnute (1987). Model theory, implementation and
testing are detailed in a manuscript currently under review.

## Installation

This package involves C++ code and the modelling framework of the
package `TMB`. Please ensure that you are able to compile source code
and that this package is running smoothly on your device. Package has
not been tested on the most recent version of R (&gt;=4.0) as of yet

And sbar package can be installed from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("lbatts/sbar")
```

## Example

A quick example of running CSA on the black-bellied anglerfish stock.
For a much more detailed look at sbar assessments run
`vignette("sbar")`. Data is structured as
[FLR](https://flr-project.org/) classes and so the appropriate syntax is
required.

``` r
library(sbar)
library(FLCore)
```

``` r
data("ank78")
data("ank78.indices")
years<-as.character(2003:2020) 
no.years<-length(years)
```

Observations needed for a CSA assessment are catch numbers and a matrix
of survey indices (catch numbers per unit effort). A survey split into a
recruit index and post-recruit index (i.e. one survey, two indices) is
the minimum requirement. In this example we use the combined IBTS survey
data that has processed already into number at age. This gives a simple
way to define recruits (age 0) and post-recruits (age 1+).

``` r
catch.no<-c(colSums(catch.n(ank78)[,years]))

no.ind = 2 
IBTS_PR_ages<-as.character(range(ank78.indices$FR_IE_IBTS)["min"]+1:range(ank78.indices$FR_IE_IBTS)["max"])
IBTS_PR_ages
#> [1] "1" "2" "3" "4" "5" "6" "7" "8" "9"

obs<-matrix(NA,nrow=no.ind,ncol=no.years)
obs[1,]<-c(colSums(index(ank78.indices$FR_IE_IBTS)["0",years],na.rm=T)) 
obs[2,]<-c(colSums(index(ank78.indices$FR_IE_IBTS)[IBTS_PR_ages,years],na.rm=T)) 

obs[obs==0]<-NA
```

Many of the settings for running CSA have defaults but the function
requires some user defined values. `indices_att` corresponds to `obs`,
indicating if indices (each row) are from the same survey (i.e. same
number) and what type of indices they are:

1.  recruit index
2.  post-recruit index
3.  undivided index

``` r
indices_att<-data.frame(survey=c(1,1),type=c(1,2))
timing=c(0.875) # survey timing
nm<-mean(m(ank78)) #natural mortality
```
