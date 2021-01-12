






schnute_datprep<-function(catch,indices,ts, mwts,tsp){

data_tmb <- list(
  obs_catch = catch,
  obs_ind = indices,
  indices_class = c(2),
  indices_ts = ts,
  mean_wts=mwts,
  nu=tsp,
  SRcode=2,
  spawn_prop=NULL)

return(data_tmb)

}

schnute_parprep_v2 <- function( q, indexsigma, B0, sigma, rho, W, f_calc, catchsigma){


params <- list(
  logq = log(q),
  logindex_sigma = log(indexsigma),
  logB0 = log(B0),
  logrec_param = log(c(1,1)),
  logitsigma = stats::qlogis(sigma),
  logrho=log(rho),
  logW=log(W),
  logf_calc = log(f_calc),
  logcatch_sigma = log(catchsigma))

return(params)

}



schnute_parprep_v1 <- function( q, indexsigma, sigma, rho, W, f_calc, catchsigma){


  params <- list(
    logq = log(q),
    logindex_sigma = log(indexsigma),
    logrec_param = log(c(1,1)),
    logitsigma = stats::qlogis(sigma),
    logrho=log(rho),
    logW=log(W),
    logf_calc = log(f_calc),
    logcatch_sigma = log(catchsigma))

  return(params)

}


