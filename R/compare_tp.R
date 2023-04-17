#' Feature expression between time-points
#'
#' This function generates a plot comparing the expression of a feature at two time-points.
#' If the size of the two data differ, bear in mind that the function will remove randomly
#' some cells to equalize the dimensions, if you want to set yourself the sampling frequencies by class
#' set "group" to true
#'
#'
#' @param Data_1 Data.frame at time-point x
#' @param Data_2 Data.frame at time-point y
#' @param Feature The feature which expression will be plotted
#' @param n The number of random selected cells (Default=100)
#' @param group Do you want to set yourself the sampling frequency of each class?
#' @param Freq Frequency of each class, in data.frame format (Var1=classes,Freq=Number of cells)
#' @param Regression Do you want to plot also the regression line?
#' @param main Title of plot
#' @param x Name of x-axis
#' @param y Name of y-axis
#'
#' @import tidymodels
#' @import reshape2
#' @import dplyr
#'
#' @return A plot comparing the expression of a feature at two time-points.
#'
#' @export
#'

compare_tp<-function(Data_1,Data_2,Feature="",n=100,group=FALSE,Freq=NULL,Regression=T,main="",x="",y=""){
  
  if (nrow(Data_1) != nrow(Data_2)){
    if (group){
      Data_1<-sub_sampling2(Freq,Data_1)
      Data_2<-sub_sampling2(Freq,Data_2)
    }
    else{
      Data_1<-Data_1[sample(nrow(Data_1), n), ]
      Data_2<-Data_2[sample(nrow(Data_2), n), ]
    }
  }
  if (Regression){
    Data_1<-scale_data(Data_1[[Feature]])
    Data_2<-scale_data(Data_2[[Feature]])
    plot(Data_1,Data_2,
         pch=16,col=c("orange","purple"),main=main,xlab=x,ylab=y)
    abline(lm(Data_2 ~ Data_1), col = 3, lwd = 3)
    legend("topright", legend=c(x, y),pch = c(16,16),
           col=c("orange", "purple"))
    p<-recordPlot()
    return(p)
  }
  else{
    Data_1<-scale_data(Data_1[[Feature]])
    Data_2<-scale_data(Data_2[[Feature]])
    plot(Data_1,Data_2,
         pch=16,col=c("orange","purple"),main=main,xlab=x,ylab=y)
    legend("topright", legend=c(x, y),pch = c(16,16),
           col=c("orange", "purple"))
    p<-recordPlot()
    return(p)
  }
}
