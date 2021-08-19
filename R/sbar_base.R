

new_sbarclass <- function(x = matrix(), survey_names,cat,ind,years) {
  stopifnot(is.matrix(x))
  
  if(any(unique(rownames(x)=="phat1"))) {
    cat_name <- "Catch numbers"
    logN<-x[rownames(x)=="lnbhat",]
    logPR<-x[rownames(x)=="lnphat",]
    logR<-x[rownames(x)=="logrhat",]
    
  }else {
    cat_name <-"Catch biomass"
    logN<-x[rownames(x)=="lnb",]
    logPR<-x[rownames(x)=="lnpr",]
    logR<-x[rownames(x)=="lnr",]
  }
  logcatch<-x[rownames(x)=="lnc",]
  logsurv<-x[rownames(x)=="logpred_survey",]
  no.years<-length(years)
  
  if(no.years!=(dim(logsurv)[1]/length(survey_names))) {
    stop("The number of survey names you have entered does not match the number of surveys used within the assessment")
  }
  
  dat<-data.frame(param=rownames(rbind(logsurv,logcatch)),lnest=rbind(logsurv,logcatch)[,1],est=exp(rbind(logsurv,logcatch)[,1]),lnse=rbind(logsurv,logcatch)[,2],upper=NA,lower=NA)
  dat$param<-as.character(dat$param)
  
  if(dim(dat)[1]==dim(logsurv)[1]){
    dat$ver<-c(rep(survey_names,times=no.years))
    dat$year<-as.numeric(c(rep(years,each=length(survey_names))))
    dat[dat$year==years[1],2:6]<-NA
    dat$obs<-c( c(ind))
  }else{
    dat$ver<-c(rep(survey_names,times=no.years),rep(cat_name,times=no.years))
    dat$year<-as.numeric(c(rep(years,each=length(survey_names)),rep(years,times=1)))
    dat$obs<-c( c(ind),cat)
  }
  dat$upper<-exp(dat$lnest + 2*dat$lnse)
  dat$lower<-exp(dat$lnest - 2*dat$lnse)
  #this should work for csa and schnute version 2 of of both
  
  
  #think about if statments here as well
  
  logF<-x[rownames(x)=="logf_calc",]
  
  dat2<-data.frame(param=rownames(rbind(logR,logPR,logN,logF)),lnest=rbind(logR,logPR,logN,logF)[,1],est=exp(rbind(logR,logPR,logN,logF)[,1]),lnse=rbind(logR,logPR,logN,logF)[,2],upper=NA,lower=NA)
  
  dat2$param<-(as.character(dat2$param))
  dat2$upper<-exp(dat2$lnest + 2*dat2$lnse)
  dat2$lower<-exp(dat2$lnest - 2*dat2$lnse)
  dat2[dat2$est==0,2:6]<-NA
  
  if(dim(dat)[1]==dim(logsurv)[1]){
    dat2$year<-as.numeric(rep(years,times=3))
  }else dat2$year<-as.numeric(rep(years,times=4))
  
  if(cat_name == "Catch numbers"){
    dat2$param[dat2$param=="lnbhat"]<-("Population numbers")
    dat2$param[dat2$param=="logf_calc"]<-("Fishing mortality")
    dat2$param[ dat2$param=="lnphat"]<-("Post-recruit numbers")
    dat2$param[ dat2$param=="logrhat"]<-("Recruitment numbers")
  }else{
    dat2$param[dat2$param=="lnb" ]<-("Population biomass")
    dat2$param[dat2$param=="logf_calc"]<-("Fishing mortality")
    dat2$param[dat2$param=="lnpr" ]<-("Post-recruit biomass")
    dat2$param[dat2$param=="lnr"]<-("Recruitment biomass")
  }
  
  dat_ls<-list(dat,dat2)
  
  structure(dat_ls, class = "sbarclass")
}
