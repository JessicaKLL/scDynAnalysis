#' Plot time step dataframe
#' 
#' Thus function generates a plot showing the dynamics of the selected feature 
#' 
#' @param data Your time step dataframe
#' @param Feature Your seleted feature
#' @param dyn_col The color of the plot
#' @param title Title of the plot
#' 
#' @import ggplot2
#' 
#' @return A plot showing the dynamics of the feature
#' 
#' @export
#' 

plot_time_series<-function(tk_df,Feature,dyn_col="red",title=""){
  colnames(data)<-c("index","value")
  p<-data %>%
    ggplot(aes(index, value)) +
    geom_line(group=1,color = dyn_col, alpha = 0.4) +
    geom_point(color = dyn_col,alpha=0.4) + 
    ylab(paste0(Feature))+xlab("Time-steps")+
    theme_classic()+theme(axis.text.x = element_blank())+
    labs(title = paste0(title))
  return(p)
}