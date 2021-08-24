#' Prepare list for Schnute adapted observation error assessment for an optimiser
#'
#' Create list with TMB framework, including data, gradients and NLL function for a Schnute adapted observation error assessment that can be optimised.
#'
#' @param version numeric, either 1, 2 or 3. This controls what deterministic equations in the model are used to derive population biomass. 1 and 2 use the fraction of of total biomass in a given year due to newly recruited fish. This fraction is derived from mean weights and detailed in the **schnute** vignette. `version = 3` is the more classical population dynamics.
#' \tabular{cc}{
#' 1 \tab whole biomass derived from recruit biomass \cr
#' 2 \tab whole biomass derived from previously exploited biomass \cr
#' 3 \tab whole biomass is a combination of recruit biomass and previously exploited biomass}
#' @param catch_b numeric vector of catch biomass over time period of assessment
#' @param indices_b matrix of biomass surveys (CPUE) of dimensions: no. of surveys x no.years
#' @param ts numeric. Survey timing parameters
#' @param mwts matrix of mean weights from sampling with dimensions: 3 x no. years. recruit mean weights \eqn{\bar{Y}} (first row), previously exploited biomass mean weights \eqn{\bar{Z}} (second row) and entire assessed biomass mean weight \eqn{\bar{X}} (third row).
#' @param tsp numeric. Timing of spawning. Default to 0 (start of year).
#' @param rho numeric. Growth parameter, slope of linear growth model.
#' @param W numeric. Growth parameter, intercept of linear growth model.
#' @param ind_l_wt numeric. Survey weighting in the likelihood. Defaults to 1 fro each survey, ie.e. equal weighting
#' @param start_q Starting values for survey catchability parameters. Default is 1e-6
#' @param start_indexsigma Starting values for survey sd parameters. Default is 0.1
#' @param start_B0 Starting parameter value for biomass at first time step. Default is 5*max(catch_b)
#' @param start_sigma Starting parameter value fraction of population that survives natural moratlity. Default is \eqn{e^{0.2}}
#' @param start_f_calc Starting parameter values for estimated fishing mortality. Default is 0.3.
#' @param start_rec_a Starting parameter value for the 'a' parameter of the Beverton-Holt stock-recruit function. The asymptotic biomass of recruits. Default is 1/5*max(catch_b).
#' @param start_rec_b Starting parameter value for the 'b' parameter of the Beverton-Holt stock-recruit function. The spawning stock biomass needed to produce a/2 on average. Default is 4*max(catch_b).
#' @param spawn_prop proportion of biomass that is mature. Defaults to 1 for each year. 
#' @param start_catchsigma Starting parameter value for catch sd. Default is 0.1
#' @param fix_sigma logical
#' @param fix_B0 logical
#' @param fix_indexsigma logical
#' @param fix_catchsigma logical
#' @param adrep logical. Whether the user would like the ADreport variables (and their derivatives) reported for starting parameters.
#' @return List with components for optimiser in R. This output is that of the function \link[TMB]{MakeADFun} from TMB
#' @export
#' @examples
# #' 


schnute_obserror<-function(version = 2,catch_b,indices_b,ts, mwts,tsp = 0, rho, W, ind_l_wt = 1, start_q = 1e-8, start_indexsigma = 0.1 ,start_B0, start_sigma = exp(-0.2) , start_f_calc = 0.3,  start_rec_a, start_rec_b, spawn_prop = 1, start_catchsigma = 0.1, fix_sigma = TRUE, fix_B0 = FALSE, fix_indexsigma = FALSE, fix_catchsigma = TRUE,adrep = FALSE){

 
  ny <- length(catch_b)
  no.survey <- dim(indices_b)[1]


  if(missing(version)) message("Argument 'version' missing. Default model version is 2")
  
  if(no.survey != length(ts)){
    stop("no. of surveys does not match ts length")
  }
  

  if(missing(ind_l_wt)){
    ind_l_wt <- rep(ind_l_wt,no.survey)
    message("Argument 'ind_l_wt' missing. Default indices likelihood weighting of 1 used for each survey")
  }
  
  if(length(ind_l_wt) == 1 & no.survey > 1){
    ind_l_wt <- rep(ind_l_wt,no.survey)
    stop("Survey weighting and number of surveys does not match")
  }else if((sum(ind_l_wt)/no.survey) != 1){
    stop("Survey weighting values must sum to 1")
  }
  
  
  if(missing(start_sigma)) {message("Argument 'start_sigma' missing. Default value used")
    start_sigma<-start_sigma}
  
  if(missing(start_B0)){ message("Argument 'start_B0' missing. Default value used")
    start_B0 <-  5*max(catch_b)
  }else start_B0 <- start_B0
  
  
  
  if(length(spawn_prop) == 1){
    spawn_prop <- rep(spawn_prop,times =ny)
    message("Argument 'spawn_prop' has length 1 . Given value used for each year")
  }else if(length(spawn_prop) == ny ){
    spawn_prop <- spawn_prop
  }else stop("spawn_prop should be length of 1 or length of number of years")
  
  

  rec_miss<-c(missing(start_rec_a),missing(start_rec_b))

  if(rec_miss[1]==TRUE) {rec_a <- (1/5)*max(catch_b)
  message("Argument 'start_rec_a' missing. Default value used")
  }else rec_a <-start_rec_a
  
  if(rec_miss[2]==TRUE) {rec_b <- 4*max(catch_b)
  message("Argument 'start_rec_b' missing. Default value used")
  }else rec_b <-start_rec_b



  if(version==2 & any(rec_miss==FALSE)) {
    warning("Argument 'version' == 2 cannot estimate recruitment parameters, starting parameters for recruitment ignored")
  }else if(version !=2 & any(rec_miss ==TRUE)){

    message("Default starting values for rec_a and/or rec_b set internally, consider setting your own")
  }






  dat_tmb<-schnute_datprep(catch_b,indices_b,ts,mwts,tsp,version,ind_l_wt, spawn_prop)
  
  if(missing(start_q)){
    start_q <- rep(start_q,no.survey)
    message("Argument 'start_q' missing. Default start q used for each survey")
  }
    
    if(length(start_q) == 1 & no.survey > 1){
      start_q <- rep(start_q,no.survey)
      message("Argument 'start q' vector is not the same length as the number of surveys. The given value will be used for each survey")
    }else if(length(start_q) > 1 & length(start_q) != no.survey) {
      stop("Error: start_indexsigma vector is not the same length as the number of surveys")
    }else start_q <-start_q
  
  
  if(missing(start_indexsigma)){
    start_indexsigma <- rep(start_indexsigma,no.survey)
    message("Argument 'start_indexsigma' missing. Default start_indexsigma used for each survey")
  }else start_indexsigma <-start_indexsigma
  
  if(length(start_indexsigma) == 1 & no.survey > 1){
    start_indexsigma <- rep(start_indexsigma,no.survey)
    message("Argument 'start_indexsigma' vector is not the same length as the number of surveys. The given value will be used for each survey")
  }else if(length(start_indexsigma) > 1 & length(start_indexsigma) != no.survey) {
    stop("Error: start_indexsigma vector is not the same length as the number of surveys")
  }else start_indexsigma <-start_indexsigma
  
  if(missing(start_catchsigma)){
    start_catchsigma <- start_catchsigma
    message("Argument 'start_catchsigma' missing. Default start_catchsigma used for each survey")
  }else start_catchsigma <-start_catchsigma
  
  if(length(start_catchsigma) > 1 ) {
    stop("Error: start_catchsigma vector has more than one value")
  }else start_catchsigma <-start_catchsigma
  
  
  if(missing(start_f_calc)) message("Argument 'start_f_calc' missing. Default value used for each year")

  if(length(start_f_calc) == 1){
    start_f_calc <- rep(start_f_calc,ny)
    message("Argument 'start_f_calc' has only one value which will be used for each year")
  }else if(length(start_f_calc) == ny){
    start_f_calc <- start_f_calc
  }else stop("Error: Mismatch in starting f_calc and number of years. Check your starting fishing mortality values")



  maprec <-ifelse(rep(version,2,)==2,c(NA,NA),c(1,2))





  par_tmb<-schnute_parprep_v2(q = start_q, indexsigma = start_indexsigma, B0 = start_B0 , sigma = start_sigma, rho, W, rec_a, rec_b, f_calc = start_f_calc, catchsigma = start_catchsigma)

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
      ADreport = adrep,
      DLL = "sbar_TMBExports")


  return(obj)



}
