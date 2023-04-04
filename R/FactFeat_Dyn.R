#' Factor-Feature dynamics
#'
#' This function generates a plot representing the dynamic change of feature and
#' factor expressions among time-points.
#'
#' @param data Entire data.frame
#' @param Features The selected features
#'
#' @import ggplot2
#' @import reshape2
#' @import tibble
#'
#' @return A plot representing the dynamic change of feature and
#' factor expressions among time-points.
#'
#' @export
#'

Fact_Feat_Dyn<-function(data,Features){
  data<-split(data,data$time_point)
  data_au<-list()
  for (i in 1:length(data)) {
    data_au[[i]]<-data[[i]][,1:length(Features)]
    data_au[[i]]$Factor<-data[[i]][,ncol(data[[i]])]
    data_au[[i]]<-colMeans(data_au[[i]])
  }
  df<-cbind(data_au[[1]],data_au[[2]])
  for (i in 3:length(data_au)) {
    df<-cbind(df,data_au[[i]])
  }
  colnames(df)<-names(data)
  df<-data.frame(df)
  df<-tibble::rownames_to_column(df,"Features")
  df<-melt(df,id.vars="Features")
  plot<-ggplot(df,aes(x=variable, y=value))+geom_smooth(aes(color=Features, group=Features))+
    theme_minimal()+theme(legend.key = element_blank())+labs(x="Time-point",y="Expression",
                                                             title="Dynamics feature-factor")
  return(plot)
}

