% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csa.R
\name{csa}
\alias{csa}
\title{Prepare list for CSA assessment to be used in optimiser}
\usage{
csa(
  catch_no,
  indices_no,
  indices_att,
  ts,
  selrec,
  start_q = 1e-08,
  start_surveycv = 0.1,
  start_prec0 = spr0,
  start_rec = srec,
  start_nmort = 0.2,
  start_f_calc = 0.3,
  start_catchcv = 0.1,
  fix_nmort = TRUE,
  fix_prec0 = FALSE,
  fix_surveycv = FALSE,
  fix_catchcv = TRUE
)
}
\arguments{
\item{catch_no}{numeric}

\item{indices_no}{matrix}

\item{indices_att}{matrix}

\item{ts}{numeric}

\item{selrec}{numeric}

\item{start_q}{numeric}

\item{start_surveycv}{numeric}

\item{start_prec0}{numeric}

\item{start_rec}{numeric}

\item{start_nmort}{numeric}

\item{start_f_calc}{numeric}

\item{start_catchcv}{numeric}

\item{fix_nmort}{logical}

\item{fix_prec0}{logical}

\item{fix_surveycv}{logical}

\item{fix_catchcv}{logical}
}
\value{
list
}
\description{
Create an object with TMB framework, including data, gradients and NLL function that can be optimised
}
