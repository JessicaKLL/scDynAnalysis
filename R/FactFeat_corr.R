#' Factor - Feature correlation
#'
#' This function calculates the correlation between each feature and the factor.
#' You can obtain the positive correlated features (with the factor) by setting
#' "positive" to true.
#'
#' @param Feat_data Data.frame of feature expressions
#' @param Fact_data Data.frame of factor expression
#' @param Features Vector of feature names
#' @param positive Obtain the positive correlated features? (Default=FALSE)
#'
#' @import corrr
#' @import tidymodels
#'
#' @return A list containing the Factor-Feature correlation or a vector
#' containing the positive correlated features.
#'
#' @export
#'

FactFeat_corr<-function(Feat_data,Fact_data,Features,positive=F){
  out_list<-list()
  for (i in 1:length(Features)) {
    out_list[[i]]<-correlate(Feat_data[,Features[i]],Fact_data)
  }
  names(out_list)<-Features
  if (positive){
    out_vec<-c()
    for (i in 1:length(Features)) {
      if(out_list[[i]][2] >= 0.3){
        out_vec<-append(out_vec,Features[i])
      }
    }
    return(out_vec)
  }
  else{
    return(out_list)
  }
}
