#' Select altered features
#'
#' This function selects the features which present a significant variance among time-points
#' by selecting features which have p-value < alpha (0.05)
#'
#' @param  data A list containing the p-values of each feature between different time-points,
#' this can be calculated the function Between_variance()
#' @param  Features Vector containing the names of the features
#'
#' @return A vector of selected features
#'
#' @export

select_alter<-function(data,Features){
  Mean_pVal<-list()
  for (i in 1:length(Features)) {
    x<-mean(data[[i]]$p.value)
    Mean_pVal[[i]]<-x
  }
  names(Mean_pVal)<-Features
  vec_select<-c()
  for (i in 1:length(Features)) {
    if(Mean_pVal[[i]] < 0.05){
      vec_select<-append(vec_select,names(Mean_pVal)[i])
    }
  }
  return(vec_select)
}
