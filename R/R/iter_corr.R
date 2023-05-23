#' Iteration of random-subsets correlation matrix
#'
#'
#' This function creates n random-subset correlation matrices between two data.frames
#'
#'
#' @param data1 Time-point data.frame 1
#' @param data2 Time-point data.frame 2
#' @param n Number of iterations (Default=100)
#' @param decreasing Show correlations decreasing? (Default=TRUE)
#' @param perc The fraction of cells selected in each cell type
#'
#' @return A list which contains n data.frames of the correlation of each feature in different time-points
#'
#' @export

iter_corr<-function(data1,data2,n=100,decreasing=TRUE,perc=0.1){
  output<-list()
  for (i in c(1:n)){
    correlation_i<-calc_correlation(data1,data2,perc=perc)
    df<-data.frame(data_1=rownames(correlation_i$correlation_matrix)[row(correlation_i$correlation_matrix)],
                   data_2=colnames(correlation_i$correlation_matrix)[col(correlation_i$correlation_matrix)],
                   corr=c(correlation_i$correlation_matrix))
    df <- df[df$data_1==df$data_2, ]
    df$data_1<-paste(df$data_1,"1",sep = "_")
    df$data_2<-paste(df$data_2,"2",sep = "_")
    df<-df[order(df$corr, decreasing = F),]
    output[[i]]<-df
  }
  return(output)
}
