#' Feature expression between time-points
#'
#' This function generates a plot comparing the expression of a feature at different time-points.
#'
#'
#' @param Data Data.frame of feature expressions
#' @param Cell_Type Cell-type of each input cell
#' @param Time_Step Time-point of each input cell.
#' @param Feature The feature which expression will be plotted
#' @param main Title of plot
#'
#' @import ggplot2
#' 
#'
#' @return A plot comparing the expression of a feature at different time-points.
#'
#' @export
#'

compare_tp<-function(Data,Cell_Type,Time_Step,Feature="",main=""){
  x<-data.frame(Data[,Feature])
  x<-data.frame(x$Data...Feature.)
  colnames(x)<-"Expression"
  x$cell_type<-Cell_Type
  x$time_point<-Time_Step
  p<-ggplot(x,aes(x=cell_type,y=Expression,col=time_point))+geom_boxplot()+
    theme_bw()+theme(axis.text.x = element_text(angle = 90))+ggtitle(main)+
    xlab("Cell Type")
  return(p)
}