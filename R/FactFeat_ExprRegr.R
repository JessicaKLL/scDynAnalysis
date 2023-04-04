#' Factor-Feature expression regression
#'
#' This function generates a list of plots comparing each feature with the same factor
#'
#' @param Feat_data Data.frame of feature expressions
#' @param Fact_data Data.frame of factor expression
#' @param Features Vector of feature names
#' @param x Name x-axis
#' @param y Name y-axis
#'
#'
#' @return A list of plots
#'
#' @export
#'

FactFeat_ExprPlot<-function(Feat_data, Fact_data, Features,y="Factor"){
  output_list<-list()
  for (i in 1:length(Features)) {
    plot(Feat_data[,Features[i]],Fact_data,
         pch=16,col=c("orange","purple"),main=paste0("Factor - ",Features[i]),
         xlab=paste0(Features[i]),
         ylab=y)
    abline(lm(Fact_data ~ Feat_data[,Features[i]]), col = 3, lwd = 3)
    legend("topright", legend=c(paste0(Features[i]), y),pch = c(16,16),
           col=c("orange", "purple"))
    p<-recordPlot()
    output_list[[i]]<-p
    plot.new()
  }
  names(output_list)<-Features
  return(output_list)
}
