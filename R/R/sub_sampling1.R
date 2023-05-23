#' Get frequency
#'
#' This function calculates the number of cells which will be randomly selected from each cell type
#' in the next step
#'
#' @param input Input data
#' @param perc The fraction of cells selected in each cell type
#'
#' @import dplyr
#'
#'
#' @return A data.frame that contains the number of cells from each cell type.
#'
#' @export
#'


sub_sampling1<-function(input,perc=0.1){
  Input<-input %>%
    group_by(cell_type) %>%
    sample_frac(perc)
  Freq<-as.data.frame(table(Input$cell_type))
  return(Freq)
}
