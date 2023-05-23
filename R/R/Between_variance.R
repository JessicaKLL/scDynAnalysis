#' Significance of the between time-point variance
#'
#' This function calculates the significance of the variance between time-points
#'
#' @param data Longitudinal data
#' @param Features Vector containing the names of the features
#' @param Time_Points Vector containing the time-point of each cell
#'
#' @return A list containing the p-values of each feature between different time-points
#'
#' @export

Between_variance<-function(data, Features, Time_Points){
  n<-length(Features)
  list_feat<-list()
  for (i in 1:n) {
    x<-Time_Points
    y<-c(data[,i])
    list_feat[[i]]<-data.frame(x,y)
    colnames(list_feat[[i]])<-c("Time_Point","data")
  }
  names(list_feat)<-Features
  Wilcox_Time<-list()
  for (i in 1:length(list_feat)) {
    Wilcox_Time[[i]]<-pairwise.wilcox.test(list_feat[[i]]$data,list_feat[[i]]$Time_Point,
                                           p.adjust.method = "holm")
    Wilcox_Time[[i]]<-broom::tidy(Wilcox_Time[[i]])
  }
  names(Wilcox_Time)<-Features
  
  return(Wilcox_Time)
}