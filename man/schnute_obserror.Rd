% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/schnute_obserror.R
\name{schnute_obserror}
\alias{schnute_obserror}
\title{Prepare an object for Schnute Adapted Observation Error Model assessment for
an optimiser}
\usage{
schnute_obserror(
  version = 2,
  catch_b,
  indices_b,
  ts,
  mwts,
  tsp = 0,
  rho,
  W,
  ind_l_wt = 1,
  start_q = 1e-08,
  start_indexsigma = 0.1,
  start_B0,
  start_sigma = exp(-0.2),
  start_f_calc = 0.3,
  start_rec_a,
  start_rec_b,
  spawn_prop = 1,
  start_catchsigma = 0.1,
  fix_sigma = TRUE,
  fix_B0 = FALSE,
  fix_indexsigma = FALSE,
  fix_catchsigma = TRUE,
  adrep = FALSE
)
}
\arguments{
\item{version}{numeric, either 1, 2 or 3. This controls what deterministic equations in the model are used to derive population biomass. 1 and 2 use the fraction of of total biomass in a given year due to newly recruited fish. This fraction is derived from mean weights and detailed in the \strong{schnute} vignette. \code{version = 3} is the more classical population dynamics.
\tabular{cc}{
1 \tab whole biomass derived from recruit biomass \cr
2 \tab whole biomass derived from previously exploited biomass \cr
3 \tab whole biomass is a combination of recruit biomass and previously exploited biomass}}

\item{catch_b}{numeric vector of catch biomass over time period of assessment}

\item{indices_b}{matrix of biomass surveys (CPUE) of dimensions: no. of surveys x no.years}

\item{ts}{numeric. Survey timing parameters}

\item{mwts}{matrix of mean weights from sampling with dimensions: 3 x no. years. recruit mean weights \eqn{\bar{Y}} (first row), previously exploited biomass mean weights \eqn{\bar{Z}} (second row) and entire assessed biomass mean weight \eqn{\bar{X}} (third row).}

\item{tsp}{numeric. Timing of spawning. Default to 0 (start of year).}

\item{rho}{numeric. Growth parameter, slope of linear growth model.}

\item{W}{numeric. Growth parameter, intercept of linear growth model.}

\item{ind_l_wt}{numeric. Survey weighting in the likelihood. Defaults to 1 fro each survey, ie.e. equal weighting}

\item{start_q}{Starting values for survey catchability parameters. Default is 1e-6}

\item{start_indexsigma}{Starting values for survey sd parameters. Default is 0.1}

\item{start_B0}{Starting parameter value for biomass at first time step. Default is 5*max(catch_b)}

\item{start_sigma}{Starting parameter value fraction of population that survives natural moratlity. Default is \eqn{e^{0.2}}}

\item{start_f_calc}{Starting parameter values for estimated fishing mortality. Default is 0.3.}

\item{start_rec_a}{Starting parameter value for the 'a' parameter of the Beverton-Holt stock-recruit function. The asymptotic biomass of recruits. Default is 1/5*max(catch_b).}

\item{start_rec_b}{Starting parameter value for the 'b' parameter of the Beverton-Holt stock-recruit function. The spawning stock biomass needed to produce a/2 on average. Default is 4*max(catch_b).}

\item{spawn_prop}{proportion of biomass that is mature. Defaults to 1 for each year.}

\item{start_catchsigma}{Starting parameter value for catch sd. Default is 0.1}

\item{fix_sigma}{logical. logical. Should survival be fixed in the model}

\item{fix_B0}{logical. Should biomass in the first time step be fixed in the model}

\item{fix_indexsigma}{logical. Should survey standard deviation be fixed in the model}

\item{fix_catchsigma}{logical. Should catch standard deviation be fixed in the model}

\item{adrep}{logical. Whether the user would like the ADreport variables (and their derivatives) reported for starting parameters.}
}
\value{
List with components for optimiser in R. This output is that of the
function \link[TMB]{MakeADFun} from TMB

Access (transformed) starting values of parameters to be estimated with
\code{x$par} (where x is the schnute_obserror object) to see what parameters are to be estimated.
}
\description{
Create an object with TMB framework, including data, gradients and NLL
function for a Schnute adapted observation error assessment that can be
optimised.
}
\details{
\code{schnute_obserror} is simply a wrapper function that gives the output
from \code{\link[TMB]{MakeADFun}} from TMB, i.e. an objective function with
derivatives, hessian etc. Otimisation and extraction of values from the
assessment models has been kept separate to allow flexibility with
optimisation methods as well as easy access to \code{\link[TMB]{MakeADFun}}
outputs such as the hessian, gradients etc. See
\url{http://kaskr.github.io/adcomp/_book/Introduction.html} and TMB
documentation for details. Users should refer to the \code{sbar} vignette by
running \code{vignette("intro_to_sbar", "sbar")} for details on the
\code{schnute_obserror} function.

The table below gives the outputs and description of the values that can be
extracted from this assessment model after optimisation and
\code{summary(TMB::sdreport(x))} of the \code{schnute_obserror} object.

\tabular{ll}{ Output         \tab Description             \cr logitq     \tab
logit transformed survey catchability              \cr logB0      \tab log
transformed biomass in the first time step       \cr logitsigma       \tab
logit transformed survival of natural mortality  \cr logindex_sigma \tab log
transformed survey standard deviation         \cr logcatch_sigma       \tab
log transformed catch standard deviation  \cr logf_calc       \tab log
transformed fishing mortality  \cr lnb           \tab log transformed total
biomass                        \cr lnpr           \tab log transformed
previously-exploited biomass               \cr lnr           \tab log
transformed recruit biomass                        \cr lnN           \tab log
transformed total numbers                        \cr lnPR           \tab log
transformed previously-exploited or post-recruit  numbers               \cr
lnR           \tab log transformed recruit numbers                        \cr
lnC           \tab log transformed predicted catch numbers
\cr biomass       \tab total biomass                        \cr B0    \tab
total biomass in the first time step  \cr N       \tab total numbers
\cr ssb       \tab spawning biomass                        \cr post_rec
\tab previously-exploited or post-recruit biomass                        \cr
PR       \tab previously-exploited or post-recruit numbers
\cr rec_bio       \tab recruit biomass                        \cr rec_no
\tab recruit numbers                        \cr C       \tab predicted catch
numbers                        \cr catch_pred       \tab predicted catch
biomass                        \cr logpred_survey \tab log transformed
predicted survey indices           \cr omega       \tab fraction of total
biomass in a given year due to newly recruited fish \cr mu       \tab
fraction of the catch removed before natural mortality \cr index_sigma
\tab Survey indices standard deviation\cr catch_sigma \tab catch standard
deviation\cr sigma          \tab survival of natural mortality\cr
logrec_param          \tab log transformed recruit parameters from a
Beverton-Holt SR function if estimated in versions 1 or 3 \cr rec_param
\tab the two recruit parameters fro a Beverton-Holt SR functionif estimated
in versions 1 or 3\cr qhat           \tab survey catchability
}
}
\examples{
\dontrun{ obj <- schnute_obserror(catch_b = catch_biomass, indices_b = obs,
ts = timing, mwts = mwts, rho = rho1, W = W1, start_sigma = sigma)}
}
