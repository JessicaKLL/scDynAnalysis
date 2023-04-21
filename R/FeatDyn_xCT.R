#' Feature dynamics per cell type
#'
#' This function generates a plot representing the dynamic change of the feature per cell type
#'
#' @param data Entire data.frame
#' @param Features The selected feature
#' @param tp_explicit Take into account time_points?
#'
#' @import ggplot2
#' @import tidyverse
#' @import tibble
#'
#' @return A plot representing the dynamic change of the feature
#'
#' @export
#'

FeatDyn_cellType<-function(data,Feature,tp_explicit=TRUE){
  df <- tibble::rownames_to_column(data, "cell")
  
  df<-split(df,df$time_point)
  df_au<-rbind(df[[1]],df[[2]])
  for (i in 3:length(df)) {
    df_au<-rbind(df_au,df[[i]])
  }
  df<-df_au
  
  x<-data.frame(df$cell,df[,Feature])
  colnames(x)<-c("Cell","Feat")
  x$time_point<-df$time_point
  x$cell_type<-df$cell_type
  x$cell_type<-as.factor(x$cell_type)
  
  if(tp_explicit){
    p<-ggplot(x,aes(x=Cell,y=Feat,group=cell_type,col=cell_type,fill=cell_type))+geom_smooth()+
      theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.line = element_line())+
      facet_grid(~time_point)+labs(title = paste0(Feature,"'s dynamics"),x="Cells",y=paste0(Feature))
    return(p)
  }
  else{
    p<-ggplot(x,aes(x=Cell,y=Feat,group=cell_type,col=cell_type,fill=cell_type))+geom_smooth()+
      theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.line = element_line())+
      labs(title = paste0("Dynamics ",Feature),x="Cells",y=paste0(Feature))
    return(p)
  }
}
