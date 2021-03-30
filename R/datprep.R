




schnute_datprep<-function(catchkg,indiceskg,ts, mwts,tsp,version,ind_l_wt){

data_tmb <- list(
  obs_catch = catchkg,
  obs_ind = indiceskg,
  indices_class = version,
  indices_ts = ts,
  mean_wts=mwts,
  nu=tsp,
  SRcode=2,
  spawn_prop=rep(1,length(catchkg)),
  l_calc_wt = c(ind_l_wt))

return(data_tmb)

}

schnute_parprep_v2 <- function( q, indexsigma, B0, sigma, rho, W,rec_a,rec_b, f_calc, catchsigma){


params <- list(
  logitq = stats::qlogis(q),
  logindex_sigma = log(indexsigma),
  logB0 = log(B0),
  logrec_param = log(c(rec_a,rec_b)),
  logitsigma = stats::qlogis(sigma),
  logrho=log(rho),
  logW=log(W),
  logf_calc = log(f_calc),
  logcatch_sigma = log(catchsigma))

return(params)

}



schnute_parprep_v1 <- function( q, indexsigma, sigma, rho, W,rec_a,rec_b){


  params <- list(
    logitq = stats::qlogis(q),
    logindex_sigma = log(indexsigma),
    logrec_param = log(c(rec_a,rec_b)),
    logitsigma = stats::qlogis(sigma),
    logrho=log(rho),
    logW=log(W))

  return(params)

}



csa_datprep<-function(catch_no, indices_no, indices_att, ts, selrec){


  data_tmb <- list(
    obs_catch = catch_no,
    obs_ind = indices_no,
    indices_class = indices_att,
    indices_ts = ts,
    sr =selrec)


  return(data_tmb)

}

csa_parprep <- function( q, surveycv, rec_est, prec0, nmort, f_calc, catchcv, srx){


  params <- list(
    logitqhat = stats::qlogis(q),
    log_surveycv = log(surveycv),
    log_catchcv = log(catchcv),
    logphat1 = (log(prec0)),
    logrhat = log(rec_est),
    logf_calc = log(f_calc),
    lognmort = log(nmort),
    logitsrx = stats::qlogis(srx)
  )
  return(params)

}


