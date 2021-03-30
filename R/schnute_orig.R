#' Prepare list for SChnute assessment in an optimiser
#'
#' Create an object with TMB framework, including data, gradients and NLL function that can be optimised
#'
#' @param version numeric
#' @param catchkg numeric
#' @param indiceskg matrix
#' @param ts numeric
#' @param mwts matrix
#' @param tsp numeric
#' @param mu numeric
#' @param rho numeric
#' @param W numeric
#' @param ind_l_wt numeric
#' @param start_q numeric
#' @param start_indexsigma numeric
#' @param start_sigma numeric
#' @param start_rec_a numeric
#' @param start_rec_b numeric
#' @param fix_sigma logical
#' @param fix_indexsigma logical
#' @return list
#' @export
#' @examples
# #' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])


schnute_orig<-function(version = 2,catchkg,indiceskg,ts, mwts,tsp = 0, mu = 0.5, rho, W, ind_l_wt = 1, start_q = 1e-8, start_indexsigma = 0.1 , start_sigma = exp(-0.2), start_rec_a, start_rec_b, fix_sigma = TRUE,fix_indexsigma = FALSE){

  sb0 <- 5*max(catchkg)
  ny <- length(catchkg)
  no.survey <- length(ts)


  if(missing(version)) message("Argument 'version' missing. Default model version is 2")

  if(missing(ind_l_wt)){
    ind_l_wt <- rep(ind_l_wt,no.survey)
    message("Argument 'ind_l_wt' missing. Default indices likelihood weighting of 1 used for each survey")
  }


  rec_miss<-c(missing(start_rec_a),missing(start_rec_b))

  if(rec_miss[1]==TRUE) {rec_a <- (1/5)*max(catchkg)
  }else rec_a <-start_rec_a

  if(rec_miss[2]==TRUE) {rec_b <- 4*max(catchkg)
  }else rec_b <-start_rec_b


  if(version==2 & any(rec_miss==FALSE)) {
    warning("Argument 'version' == 2 cannot estimate recruitment parameters, starting parameters for recruitment ignored")
  }else if(version !=2 & any(rec_miss ==TRUE)){

    message("Default starting values for rec_a and/or rec_b set internally, consider setting your own")
  }




  dat_tmb<-schnute_datprep(catchkg,indiceskg,ts,mwts,tsp,version,ind_l_wt)

  if(missing(start_q)){
    start_q <- rep(start_q,no.survey)
    message("arg: 'start_q' missing. Default start q used for each survey")
  }else


  if(missing(start_indexsigma)){
    start_indexsigma <- rep(start_indexsigma,no.survey)
    message("arg: 'start_indexsigma' missing. Default start_indexsigma used for each survey")
  }else start_indexsigma <-start_indexsigma




  maprec <-ifelse(rep(version,2,)==2,c(NA,NA),c(1,2))

  par_tmb<-schnute_parprep_v1(q = start_q, indexsigma = start_indexsigma, sigma = start_sigma, rho, W, rec_a, rec_b)

  dat_tmb$mu <- mu

    obj <- TMB::MakeADFun(
      data = c(model = "schnute_orig",dat_tmb),
      parameters = par_tmb,
      map = list(
        logindex_sigma = factor(ifelse(rep(fix_indexsigma,no.survey)==T,NA,c(seq(from=1, to=no.survey, by = 1)))),
        logW = factor(NA),
        logrho = factor(NA),
        logrec_param = factor(maprec),
        logitsigma = factor(ifelse(fix_sigma==T,NA,1))
      ),
      hessian = TRUE,
      silent = TRUE,
      DLL = "sbar_TMBExports")



  return(obj)



}
