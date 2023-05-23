#' Scale the data
#' 
#' Scale the dataset
#' 
#' @return Scaled data
#' 
#' @export
#' 

scale_data<-function(data){
  x_mean<-mean(data)
  x_sd<-sd(data)
  x_scale<-data
  for (i in 1:length(x_scale)) {
    x_scale[i]<-(x_scale[i]-x_mean)/x_sd
  }
  return(x_scale)
}
