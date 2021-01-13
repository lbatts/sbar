




schnute_datprep<-function(catch,indices,ts, mwts,tsp){

data_tmb <- list(
  obs_catch = catch,
  obs_ind = indices,
  indices_class = c(2),
  indices_ts = ts,
  mean_wts=mwts,
  nu=tsp,
  SRcode=2,
  spawn_prop=NA)

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



csa_datprep<-function(catch_no, indices, indices_att, ts, sr){


  data_tmb <- list(
    obs_catch = catch_no,
    obs_ind = indices,
    indices_class = indices_att,
    indices_ts = ts,
    sr =sr)


  return(data_tmb)

}

csa_parprep <- function( q, surveycv, rec, prec0, nmort, f_calc, catchcv, srx){


  params <- list(
    logitqhat = stats::qlogis(q),
    log_surveycv = log(surveycv),
    log_catchcv = log(catchcv),
    logphat1 = (log(prec0)),
    logrhat = log(rec),
    logf_calc = log(f_calc),
    lognmort = log(nmort),
    logitsrx = stats::qlogis(srx)
  )
  return(params)

}


