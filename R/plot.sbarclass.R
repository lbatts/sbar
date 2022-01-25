#' Plot object of with class 'sbarclass'
#'
#' Plot an object with class:'sbarclass', plot fit or stock predictions with ggplot functions
#' 
#' See vignette for a detailed example of this function in use.
#'
#' @param x object of class:sbarclass
#' @param out  a character that specifies whether the "fit" or "stock" should be plotted
#' @param ... No specific usage
#' @return list
#' @seealso \code{\link{makesbarclass}} for to conversion into an object that can be used by this function.
#' @examples \dontrun{ 
#' plot(x,out="fit")
#'   plot(x,out="stock")}
#' @exportS3Method 
#' 

plot.sbarclass<-function(x,out="fit",...){
 
  stopifnot(class(x)=="sbarclass")
  if(out=="fit"){
    datgg<-x[[1]]
    s=.5
    t=1
    p2<-ggplot2::ggplot(data=datgg, ggplot2::aes_string(x="year", y="est"))+ ggplot2::geom_ribbon(ggplot2::aes_string(ymin="lower", ymax="upper"),fill="lightgrey",color="lightgrey")+ggplot2::geom_line(ggplot2::aes_(colour=quote("Predicted index")),size=s) +ggplot2::geom_point(ggplot2::aes_(y = quote(obs), col= quote("Observed index")),size=t,shape=8)+ ggplot2::geom_line()+ggplot2::geom_line(ggplot2::aes_(y = quote(obs), col= quote("Observed index")),size=s,linetype="dashed")+ ggplot2::scale_colour_manual(name='', values=c('Predicted index'='black', 'Observed index'='black'))+ggplot2::scale_fill_manual(name='', values=c('Predicted index'=1, 'Observed index'=NA))+ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(2,1), shape=c(8,NA ))))+ggplot2::theme(text = ggplot2::element_text(size=20))
    
    p2<-p2+ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),strip.background =ggplot2::element_rect(color="black",fill="lightgrey"))+ggplot2::labs(y="", x="Year")+ggplot2::facet_wrap(~ver,scales="free_y",ncol = 2)

  }else{
    datgg<-x[[2]]
    p2<-ggplot2::ggplot(data=datgg, ggplot2::aes_string(x="year", y="est"))+ ggplot2::geom_ribbon(ggplot2::aes_string(ymin="lower", ymax="upper"),fill="lightgrey",color="lightgrey")+ ggplot2::geom_line()
    p2<-p2+ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),strip.background =ggplot2::element_rect(color="black",fill="lightgrey"))+ggplot2::labs(y="", x="Year")+ggplot2::facet_wrap(~param,scale="free_y",ncol = 2)
  }

  p2
}
