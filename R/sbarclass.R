#' Convert sdreport summary into an object of for plotting with class: sbarclass
#'
#' Create an object with class:'sbarclass', which can then be used to quickly plot fit and stock predictions from an sbar assessment
#' 
#' See vignette for a detailed example of this function in use.
#' 
#' @param x the matrix output from \code{summary(sdreport(x))} , where \code{x} is an optimised sbar stock assessment model.
#' @param survey_names character vector of survey names
#' @param cat numeric vector of catch that was used in the sbar assessment. Used to plot the fit of the assessment
#' @param ind matrix of survey indices that were used in the sbar assessment. Used to plot the fit of the assessment
#' @param years numeric vector of years
#' @return An object of class:sbarclass
#' @seealso \code{\link{plot.sbarclass}} to plot a sbar assessment after using this function.
#' @export
#' @examples \dontrun{ survnames<- c("IBTS recruits (CPUE)","IBTS post-recruits (CPUE)")
#' x <- makesbarclass(obs.srep,survnames,catch.no,obs,years)}

makesbarclass <- function(x = matrix(), survey_names,cat,ind,years) {
  x <- as.matrix(x)
  survey_names<-as.character(survey_names)
  ind<-as.matrix(ind)
  years<-as.numeric(years)
  
  new_sbarclass(x, survey_names,cat,ind,years)
}

