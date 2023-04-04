#' Sub-sampling
#'
#' This function select randomly n cells from each cell_type according to the input frequencies
#'
#'
#' @param data The data which will be sub_sampled
#' @param Freq The number of samples that are selected from each group
#'
#' @import dplyr
#' @import tidymodels
#' @import tidyr
#' @import purrr
#'
#'
#' @return A sub-sampled input data.frame
#'
#' @export
#'

sub_sampling2<-function(Freq,data){
  nested_data<- data %>%
    rownames_to_column('cell') %>%
    group_by(cell_type) %>%
    nest() %>%
    ungroup() %>%
    mutate(n=c(Freq$Freq))
  sampled_data<-nested_data %>%
    mutate(samp=map2(data,n,sample_n))
  sampled_data<-sampled_data %>%
    select(-data) %>%
    unnest(samp) %>%
    column_to_rownames('cell')
  sampled_data$n<-NULL
  return(sampled_data)
}

