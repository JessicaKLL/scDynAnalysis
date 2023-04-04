#' Factor - Feature expression comparison
#'
#' This function generates several plots comparing the expression of each feature
#' and the factor taking into account the parameters time-point and cell-type.
#'
#' @param Feat_data Data.frame of feature expressions
#' @param Fact_data Data.frame of factor expression
#' @param Features Vector of feature names
#'
#' @import ggplot2
#'
#' @return A list containing the plots
#'
#' @export
#'

FactFeat_ExprPlot<-function(Feat_data,Fact_data,Features){
  output_list<-list()
  for (i in 1:length(Features)) {
    x<-Feat_data[,Features[i]]
    y<-Fact_data
    cell_type<-Feat_data$cell_type
    time_point<-Feat_data$time_point
    data<-data.frame(x,y,cell_type,time_point)
    plot<-ggplot(data)+geom_jitter(aes(x=x,y=y,col=cell_type))+
      facet_grid(cols=vars(time_point))+
      labs(title = paste0("Factor - ",Features[i]),x=paste0(Features[i]),y="Factor")
    output_list[[i]]<-plot
  }
  names(output_list)<-Features
  return(output_list)
}
