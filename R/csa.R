#' Prepare list for CSA assessment to be used in optimiser
#'
#' Create an object with TMB framework, including data, gradients and NLL function that can be optimised
#'
#' @param catch_no numeric
#' @param indices matrix
#' @param indices_att matrix
#' @param ts numeric
#' @param sr numeric
#' @param start_q numeric
#' @param start_surveycv numeric
#' @param start_prec0 numeric
#' @param start_nmort numeric
#' @param start_f_calc numeric
#' @param start_catchcv numeric
#' @param fix_nmort logical
#' @param fix_prec0 logical
#' @param fix_surveycv logical
#' @param fix_catchcv logical
#' @return list
#' @export
#' @examples
# #' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])


csa<-function(catch_no,indices,indices_att,ts, sr, start_q = 1e-8, start_surveycv = 0.1 ,start_prec0 = spr0, start_rec = srec, start_nmort = 0.2 , start_f_calc = 0.3, start_catchcv = 0.1, fix_nmort = TRUE, fix_prec0 = FALSE, fix_surveycv = FALSE, fix_catchcv = TRUE){

  spr0 <- 4*max(catch_no)
  srec <- 2*max(catch_no)
  ny <- length(catch_no)
  no.index <- dim(indices)[1]
  no.survey <- length(unique(indices_att[,1]))

  tmp<-as.data.frame(indices_att)
  colnames(tmp) <- c("survey","type")
  tmp<-reshape2::dcast(indices_att,survey~type,value.var = "type")
  tmp$test <- ifelse(rowSums(tmp[2:3])==3,0,NA)
  (sum(tmp$test,na.rm = T)==0)


  if((sum(tmp$test,na.rm = T)!=0)){
    stop("Assessment needs minimum of one survey with both a recruit (class 1) and post-recruit (class 2) indices")
  }


  if(no.survey != length(ts)){
    stop("no. of surveys does not match ts length")
  }

  message(paste("Summary: Input is",no.survey,"unique survey(s)","but",no.index,"indices"))

indices_att <- as.matrix(indices_att)

  dat_tmb<-csa_datprep(catch_no,indices, indices_att, ts, sr)

  if(length(start_q) == 1 & no.survey > 1){
    start_q <- rep(start_q,no.survey)
    message("default start q used for each survey")
  }else start_q <-start_q


  if(length(start_surveycv) == 1 & no.survey > 1){
    start_surveycv <- rep(start_surveycv,no.survey)
    message("default start_surveycv used for each survey")
  }else start_surveycv <-start_surveycv

  if(length(start_f_calc) == 1){
    start_f_calc <- rep(start_f_calc,ny)
  }else if(length(start_f_calc) == ny){
    start_f_calc <- start_f_calc
  }else stop("start_f_calc should be length of 1 or number of years")

  if(length(start_rec) == 1){
    start_rec <- rep(start_rec,ny)
  }else if(length(start_rec) == ny){
    start_rec <- start_rec
  }else stop("start_rec should be length of 1 or number of years")



    par_tmb<-csa_parprep(q = start_q, surveycv = start_surveycv, prec0 = start_prec0 , rec = start_rec, nmort = start_nmort, f_calc = start_f_calc, catchcv = start_catchcv,srx = rep(1,ny))


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
      DLL = "sbar_TMBExports")


  return(obj)



}
