#' Convert sdreport summary into an object of for plotting with class: sbarclass
#' #'
#' Create an object with class:'sbarclass', which can then be used to quickly plot fit and stock predictions from an sbar assessment
#'
#' @param x matrix
#' @param survey_names character
#' @param cat numeric
#' @param ind matrix
#' @param years numeric
#' @return list
#' @export
#' @examples
# #' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])

makesbarclass <- function(x = matrix(), survey_names,cat,ind,years) {
  x <- as.matrix(x)
  survey_names<-as.character(survey_names)
  ind<-as.matrix(ind)
  years<-as.numeric(years)
  
  new_sbarclass(x, survey_names,cat,ind,years)
}

