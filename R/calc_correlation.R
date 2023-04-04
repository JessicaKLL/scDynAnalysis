#' Calculate the correlation
#'
#' This function calculates the correlation between the sub-sampled two data.frames
#'
#' @param data1 Time-point data.frame 1
#' @param data2 Time-point data.frame 2
#' @param perc The fraction of cells selected in each cell type
#'
#' @import dplyr
#' @import tidymodels
#'
#' @return A list which contains the cells selected in data1 and data2, and a correlation matrix.
#'
#' @export

calc_correlation<-function(data1, data2,perc=0.1){
  Filter<-sub_sampling1(data1,perc = perc)
  data_1<-sub_sampling2(Filter,data1)
  data_2<-sub_sampling2(Filter,data2)
  list_of_results<-list(data1=rownames(data_1),
                        data2=rownames(data_2))
  data_1$cell_type<-NULL
  data_1$n<-NULL
  data_2$cell_type<-NULL
  data_2$n<-NULL
  corr_m<-cor(data_1,data_2)
  list_of_results_updated <- c(list_of_results, correlation_matrix = list(corr_m))
  return(list_of_results_updated)
}
