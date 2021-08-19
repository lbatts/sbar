#' Plot object of with class 'sbarclass'
#'
#' Plot object with class:'sbarclass', plot fit or stock predictions 
#'
#' @param x sbar
#' @param out character
#' @param ... No specific usage
#' @return list
#' @examples
#' @exportS3Method 
# #' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])

plot.sbarclass<-function(x,out="fit",...){
 
  stopifnot(class(x)=="sbarclass")
  if(out=="fit"){
    datgg<-x[[1]]
    s=.5
    t=1
    p2<-ggplot2::ggplot(data=datgg, ggplot2::aes(x=year, y=est))+ ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper),fill="lightgrey",color="lightgrey")+ggplot2::geom_line(ggplot2::aes(colour="Predicted index"),size=s) +ggplot2::geom_point(ggplot2::aes(y = obs, col= "Observed index"),size=t,shape=8)+ ggplot2::geom_line()+ggplot2::geom_line(ggplot2::aes(y = obs, col= "Observed index"),size=s,linetype="dashed")+ ggplot2::scale_colour_manual(name='', values=c('Predicted index'='black', 'Observed index'='black'))+ggplot2::scale_fill_manual(name='', values=c('Predicted index'=1, 'Observed index'=NA))+ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype=c(2,1), shape=c(8,NA ))))+ggplot2::theme(text = ggplot2::element_text(size=20))
    
    p2<-p2+ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),strip.background =ggplot2::element_rect(color="black",fill="lightgrey"))+ggplot2::labs(y="", x="Year")+ggplot2::facet_wrap(~ver,scales="free_y",ncol = 2)

  }else{
    datgg<-x[[2]]
    p2<-ggplot2::ggplot(data=datgg, ggplot2::aes(x=year, y=est))+ ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper),fill="lightgrey",color="lightgrey")+ ggplot2::geom_line()
    p2<-p2+ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),strip.background =ggplot2::element_rect(color="black",fill="lightgrey"))+ggplot2::labs(y="", x="Year")+ggplot2::facet_wrap(~param,scale="free_y",ncol = 2)
  }

  p2
}
