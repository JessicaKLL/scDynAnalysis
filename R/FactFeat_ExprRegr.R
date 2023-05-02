#' Factor-Feature expression regression
#'
#' This function generates a plot comparing a feature with a factor
#'
#' @param Feat_data Data.frame of feature expressions
#' @param Fact_data Data.frame of factor expression
#' @param Cell_Type Cell-type of each input cell
#' @param Time_Step Split the plot by time-points? If 'yes', input the time-point of each input cell.
#' Otherwise, ignore this parameter.
#' @param Feature Which feature do you want to check?
#' @param main The title of the plot
#'
#' @import ggplot2
#' 
#' @return A plot comparing the expression of the feature and the factor.
#'
#' @export
#'

FactFeat_ExprRegr<-function(Feat_data,Fact_data,Cell_Type,Time_Step=NULL,Feature="",main=""){
  x<-data.frame(Feat_data[,Feature])
  colnames(x)<-"Expression"
  x$FeatFact<-Feature
  x$cell_type<-Cell_Type
  if(missing(Time_Step)){
    y<-data.frame(Fact_data)
    colnames(y)<-"Expression"
    y$FeatFact<-"FACTOR"
    y$cell_type<-Cell_Type
    z<-rbind(x,y)
    p<-ggplot(z,aes(x=cell_type,y=Expression,col=FeatFact))+geom_boxplot()+
      theme_bw()+theme(axis.text.x = element_text(angle = 90))+
      scale_color_manual(values = c("purple","orange"))+ggtitle(main)+
      xlab("Cell Type")
    return(p)
  }
  else{
    x$time_point<-Time_Step
    y<-data.frame(Fact_data)
    colnames(y)<-"Expression"
    y$FeatFact<-"FACTOR"
    y$cell_type<-Cell_Type
    y$time_point<-Time_Step
    z<-rbind(x,y)
    p<-ggplot(z,aes(x=cell_type,y=Expression,col=FeatFact))+geom_boxplot()+facet_grid(~time_point)+
      theme_bw()+theme(axis.text.x = element_text(angle = 90))+
      scale_color_manual(values = c("purple","orange"))+ggtitle(main)+
      xlab("Cell Type")
    return(p)
  }
}