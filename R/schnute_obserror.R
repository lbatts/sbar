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
#' @param rho numeric
#' @param W numeric
#' @param ind_l_wt numeric
#' @param start_q numeric
#' @param start_indexsigma numeric
#' @param start_B0 numeric
#' @param start_sigma numeric
#' @param start_f_calc numeric
#' @param start_rec_a numeric
#' @param start_rec_b numeric
#' @param start_catchsigma numeric
#' @param fix_sigma logical
#' @param fix_B0 logical
#' @param fix_indexsigma logical
#' @param fix_catchsigma logical
#' @return list
#' @export
#' @examples
# #' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])


schnute_obserror<-function(version,catchkg,indiceskg,ts, mwts,tsp = 0, rho, W, ind_l_wt = 1, start_q = 1e-8, start_indexsigma = 0.1 ,start_B0 = sb0, start_sigma = exp(-0.2) , start_f_calc = 0.3,  start_rec_a, start_rec_b, start_catchsigma = 0.1, fix_sigma = TRUE, fix_B0 = FALSE, fix_indexsigma = FALSE, fix_catchsigma = TRUE){

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
    message("Argument 'start_q' missing. Default start q used for each survey")
  }else


    if(missing(start_indexsigma)){
      start_indexsigma <- rep(start_indexsigma,no.survey)
      message("Argument 'start_indexsigma' missing. Default start_indexsigma used for each survey")
    }else start_indexsigma <-start_indexsigma

  if(missing(start_f_calc)) message("Argument 'start_f_calc' missing. Default value used for each year")

  if(length(start_f_calc) == 1){
    start_f_calc <- rep(start_f_calc,ny)
  }else if(length(start_f_calc) == ny){
    start_f_calc <- start_f_calc
  }else stop("Error: Mismatch in starting f_calc and number of years. Check your starting fishing mortality values")



  maprec <-ifelse(rep(version,2,)==2,c(NA,NA),c(1,2))





  par_tmb<-schnute_parprep_v2(q = start_q, indexsigma = start_indexsigma, B0 = sb0 , sigma = start_sigma, rho, W, rec_a, rec_b, f_calc = start_f_calc, catchsigma = start_catchsigma)

    obj <- TMB::MakeADFun(
      data = c(model = "schnute_new_V2",dat_tmb),
      parameters = par_tmb,
      map = list(
        logB0=factor(ifelse(fix_B0==T,NA,1)),
        logindex_sigma = factor(ifelse(rep(fix_indexsigma,no.survey)==T,NA,c(seq(from=1, to=no.survey, by = 1)))),
        #logf_calc= factor(rep(NA,ny)),
        logitsigma = factor(ifelse(fix_sigma==T,NA,1)),
        logcatch_sigma = factor(ifelse(fix_catchsigma==T,NA,1)),
        logrec_param = factor(maprec),
        logW = factor(NA),
        logrho = factor(NA)
      ),
      hessian = TRUE,
      silent = TRUE,
      DLL = "sbar_TMBExports")


  return(obj)



}
