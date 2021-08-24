#' Prepare list for CSA assessment to be used in optimiser
#'
#' Create list with TMB framework, including data, gradients and NLL function for a CSA assessment that can be optimised.
#'
#' @param catch_no numeric vector of catch numbers over time period of assessment
#' @param indices_no matrix of survey indices (numbers) of dimensions: no. of indices x no.years
#' @param indices_att matrix of survey indices attributes of dimensions: no. of indices x 2. First column defines survey and second column defines survey type (1 = recruit index, 2 post-recruit index, 3 = whole asessed population index). For example the minimum needed to run CSA is one survey split into a recruit index and a post-recruit index, the attribute matrix should look like: \tabular{cc}{
#' survey \tab type\cr
#' 1 \tab 1\cr
#' 1 \tab 2}

#' @param ts numeric. Survey timing parameters
#' @param selrec matrix of selection proportions of the recruit indices (in comparison to the post-recruit index) if known. Dimensions: no. of recruit indices x no. years. Defaults to 1 for all years (i.e. no difference between recruit and post-recruit indices of a survey)
#' @param start_q Starting values for survey catchability parameters. Default is 1e-6
#' @param start_surveycv Starting values for survey cv parameters. Default is 0.1
#' @param start_prec0 Starting parameter value for post-recruit numbers at first time step. Default is 4*max(catch.no).
#' @param start_rec Starting parameter values for estimated recruit numbers. Default is 2*max(catch.no).
#' @param start_nmort Starting parameter value for natural mortality. Default is 0.2
#' @param start_f_calc Starting parameter values for estimated fishing mortality. Default is 0.3.
#' @param start_catchcv  Starting parameter value for catch cv. Default is 0.1
#' @param fix_nmort logical
#' @param fix_prec0 logical
#' @param fix_surveycv logical
#' @param fix_catchcv logical
#' @return List with components for optimiser in R. This output is that of the function \link[TMB]{MakeADFun} from TMB
#' @export
#' @example 
# #' 


csa<-function(catch_no,indices_no,indices_att,ts, selrec = 1, start_q = 1e-8, start_surveycv = 0.1 ,start_prec0, start_rec, start_nmort = 0.2 , start_f_calc = 0.3, start_catchcv = 0.1, fix_nmort = TRUE, fix_prec0 = FALSE, fix_surveycv = FALSE, fix_catchcv = TRUE){

  
  
  ny <- length(catch_no)
  no.index <- dim(indices_no)[1]
  no.survey <- length(unique(indices_att[,1]))

  tmp<-as.data.frame(indices_att)
  colnames(tmp) <- c("survey","type")
  tmp<-reshape2::dcast(indices_att,survey~type,value.var = "type")
  tmp$test <- ifelse(rowSums(tmp[2:3])==3,0,NA)
  
  nr<-length(indices_att[indices_att[,2]==1,2])


  if((sum(tmp$test,na.rm = T)!=0)){
    stop("Assessment needs minimum of one survey with both a recruit (class 1) and post-recruit (class 2) indices")
  }


  if(no.survey != length(ts)){
    stop("no. of surveys does not match ts length")
  }
  
  if(missing(selrec)) {
    message("Argument 'selrec' missing. Recruits index/indices assumed fully selected")
    selrec <- matrix(selrec,nrow=nr,ncol=ny)
  }
  
  if(length(selrec) == 1){
    selrec <- matrix(selrec,nrow=nr,ncol=ny)
    message("Argument 'selrec' has length 1 . Given value used for each year")
  }else if(dim(selrec)[1] == nr & dim(selrec)[2] == ny ){
    selrec <- as.matrix(selrec)
  }else stop("selrec should be length of 1 or matrix of dim number of recruit indices x number of years")
  
  
  #message(paste("Summary: Input is",no.survey,"unique survey(s)","but",no.index,"indices"))

indices_att <- as.matrix(indices_att)

  dat_tmb<-csa_datprep(catch_no,indices_no, indices_att, ts, selrec)

  if(missing(start_q)){
    start_q <- rep(start_q,no.survey)
    message("Argument 'start_q' missing. Default start q used for each survey")
  }
  
  if(length(start_q) == 1 & no.survey > 1){
    start_q <- rep(start_q,no.survey)
    message("start q vector is not the same length as the number of surveys. The given value will be used for each survey")
  }else if(length(start_q) > 1 & length(start_q) != no.survey) {
    stop("Error: start_q vector is not the same length as the number of surveys")
  }else start_q <-start_q
  
  if(missing(start_surveycv)) {
    message("Argument 'start_surveycv' missing. Default value used for each survey")
    start_surveycv <- rep(start_surveycv,no.survey)
  }
  
  if(length(start_surveycv) == 1 & no.survey > 1){
    start_surveycv <- rep(start_surveycv,no.survey)
    message("start_surveycv vector is not the same length as the number of surveys. The given value will be used for each survey")
  }else if(length(start_surveycv) > 1 & length(start_surveycv) != no.survey) {
    stop("Error: start_surveycv vector is not the same length as the number of surveys")
  }else start_surveycv <-start_surveycv
  
  
  if(missing(start_catchcv)) {
    message("Argument 'start_catchcv' missing. Default value used for each survey")
    start_catchcv <- start_catchcv
  }
  
 if(length(start_catchcv) > 1) {
    stop("Error: start_catchcv vector has more than one value")
  }else start_catchcv <-start_catchcv
  
  
  
  if(missing(start_f_calc)) {
    message("Argument 'start_f_calc' missing. Default value used for each year")
    start_f_calc <- rep(start_f_calc,ny)
  }
  
  if(length(start_f_calc) == 1){
    start_f_calc <- rep(start_f_calc,ny)
    message("Argument 'start_f_calc' has length 1 . Given value used for each year")
  }else if(length(start_f_calc) == ny){
    start_f_calc <- start_f_calc
  }else stop("start_f_calc should be length of 1 or number of years")

  if(missing(start_prec0)){ message("Argument 'start_prec0' missing. Default value used")
    start_prec0 <- 4*max(catch_no)
  }else start_prec0 <- start_prec0 
  
  
  if(missing(start_rec)){ message("Argument 'start_rec' missing. Default value used for each year")
  start_rec <- rep((2*max(catch_no)),ny)
}
  
  if(length(start_rec) == 1){
    start_rec <- rep(start_rec,ny)
    message("Argument 'start_rec' has length 1 . Given value used for each year")
  }else if(length(start_rec) == ny){
    start_rec <- start_rec
  }else stop("start_rec should be length of 1 or number of years")

  
  if(missing(start_nmort)) {message("Argument 'start_nmort' missing. Default value used")
    start_nmort<-start_nmort}
  

    par_tmb<-csa_parprep(q = start_q, surveycv = start_surveycv, prec0 = start_prec0 , rec_est = start_rec, nmort = start_nmort, f_calc = start_f_calc, catchcv = start_catchcv,srx = rep(1,ny))


    obj <- TMB::MakeADFun(
      data = c(model = "csa",dat_tmb),
      parameters = par_tmb,
      map = list(
        logphat1=factor(ifelse(fix_prec0==T,NA,1)),
        log_surveycv = factor(ifelse(rep(fix_surveycv,no.survey)==T,NA,c(seq(from=1, to=no.survey, by = 1)))),
        log_catchcv = factor(ifelse(fix_catchcv==T,NA,1)),
        lognmort = factor(ifelse(fix_nmort==T,NA,1)),
        logitsrx = factor(rep(NA,ny))#,
        ),
      hessian = TRUE,
      silent = TRUE,
      DLL = "sbar_TMBExports")

    if(!is.finite(obj$fn(obj$par))){

      stop("Objective function is not finite so will not work in an optimiser. Try diffrent starting parameters.")
    
  }else return(obj)



}
