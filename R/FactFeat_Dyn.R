#' Factor-Feature dynamics
#'
#' This function generates a plot representing the dynamic change of feature and
#' factor expressions.
#'
#' @param data Entire data.frame
#' @param Features The selected features and the factor
#' @param Fact_data Factor expression data
#' @param tp_explicit Take into account time_points?
#'
#' @import ggplot2
#' @import tidyverse
#' @import tibble
#'
#' @return A plot representing the dynamic change of feature and
#' factor expressions.
#'
#' @export
#'

FactFeat_Dyn<-function(data,Features,Fact_data,tp_explicit=TRUE){
  df <- tibble::rownames_to_column(data, "cell")
  df<-split(df,df$time_point)
  df_au<-rbind(df[[1]],df[[2]])
  for (i in 3:length(df)) {
    df_au<-rbind(df_au,df[[i]])
  }
  df<-df_au
  x<-data.frame(df$cell,df[,Features[1]])
  colnames(x)<-c("Cell","Feature")
  x$time_point<-df$time_point
  x$Features<-Features[1]
  for (i in 2:length(Features)) {
    x_i<-data.frame(df$cell,df[,Features[i]])
    colnames(x_i)<-c("Cell","Feature")
    x_i$time_point<-df$time_point
    x_i$Features<-Features[i]
    x<-rbind(x,x_i)
  }
  x$FeatFact<-"Feature"
  x_i<-data.frame(df$cell,Fact_data)
  colnames(x_i)<-c("Cell","Feature")
  x_i$time_point<-df$time_point
  x_i$Features<-"FACTOR"
  x_i$FeatFact<-"Factor"
  Features<-append(Features,"FACTOR")
  x<-rbind(x,x_i)
  x$FeatFact

  if(tp_explicit){
    p<-ggplot(x,aes(Cell,Feature,group=Features,col=Features,fill=stage(Features, after_scale = alpha(fill, .05))))+geom_smooth()+
      theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.line = element_line(),
            panel.grid = element_blank(),panel.background = element_blank())+
      facet_grid(~time_point)+labs(title = "Feature's dynamics",x="Cells",y="Features")
    return(p)
  }
  else{
    p<-ggplot(x,aes(x=Cell,y=Feature,group=Features,col=Features,fill=stage(Features, after_scale = alpha(fill, .05))))+geom_smooth()+
      theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.line = element_line(),
            panel.grid = element_blank(),panel.background = element_blank())+
      labs(title = "Feature's dynamics",x="Cells",y="Features")
    return(p)
  }
}
