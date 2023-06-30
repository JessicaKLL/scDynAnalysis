#' Generate time step data
#' 
#' This function takes your longitudinal dataframe and generates a time step dataframe for a selected feature
#' 
#' @param data Your data
#' @param Feature Your selected feature
#' 
#' @import dplyr
#' 
#' @return A time step dataframe, where "index" are the time steps
#' 
#' @export
#' 


Gen_tk_data<-function(data,Feature){
  data <- data %>%
    tk_tbl() %>%
    mutate(index = clusters)
  data_idx<-data$index
  data<-data[,Feature]
  data<-cbind(data_idx,data)
  colnames(data)<-c("index",Feature)
  return(data)
}
