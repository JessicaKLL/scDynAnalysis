#' Significance of the within time-point variance
#'
#' This function calculates the significance of the variance within time-points
#'
#' @param data1 Time-point data.frame 1
#' @param data2 Time-point data.frame 2
#' @param n Number of iterations (Default=100)
#'
#' @import tidymodels
#'
#' @return A data.frame containing the p-values between the iterations
#'
#' @export
#' 


Within_variance<-function(data1,data2,n=100){
  Output_list<-list()

  z<-iter_corr(data1,data2,n=n,decreasing=F)

  for (i in c(1:n)) {
    id<-rep(i,nrow(z[[i]]))
    x<-data.frame(Iteration=id,Correlation=z[[i]])
    Output_list[[i]]<-x
  }

  Output_df<- Reduce(function(x, y) merge(x, y, all=TRUE), Output_list)

  Wilcoxon<-pairwise.wilcox.test(ADT_Day0_Day2_corr_df$Correlation,
                                             ADT_Day0_Day2_corr_df$Iteration, p.adjust.method="holm")
  Wilcoxon<-broom::tidy(Wilcoxon)

  return(Wilcoxon)
}

