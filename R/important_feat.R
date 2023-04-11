#' Important features
#'
#' This function selects n features by their importance calculated by the random forest model
#'
#' @param model The random forest model
#' @param data The dataframe used to generate the model
#' @param n The number of features that will be selected (Default=30)
#' @param feat_num The number of features
#' @param Decrease In which order should the features be selected? Increasing or decreasing?
#' (Default=TRUE)
#' @param main The title of the feature importance plot
#'
#' @import ggplot2
#' @import dplyr
#'
#' @return List containing the selected features and a feature importance plot
#'
#' @export

important_feat<-function(model,data,n=30,feat_num,Decrease=TRUE,main=""){
  cell_num<-length(unique(data$cell_type))
  imp_data<-data.frame(model$importance[,1:cell_num])
  features<-colnames(data[,1:feat_num])
  cell_type<-colnames(imp_data)
  Importance<-as.vector(imp_data)
  imp_df<-data.frame(features,Importance)
  colnames(imp_df)<-c("Features","Importance")
  imp_df$Importance<-as.numeric(imp_df$Importance)
  rownames(imp_df)<-NULL
  imp_df<-imp_df[order(imp_df$Importance, decreasing=Decrease),]
  imp_df<-imp_df[1:n,]
  Selected_features<-imp_df$Features
  plot<-ggplot(imp_df, aes(x=reorder(Features,Importance), y=Importance, fill=Importance))+
    geom_bar(stat = "identity",position="dodge")+coord_flip()+
    ylab("Variable Importance") + xlab("") +
    ggtitle(main)+
    scale_fill_gradientn(colors = c("yellow","orange","purple"))
  output_list<-list()
  output_list[[1]]<-Selected_features
  output_list[[2]]<-plot
  names(output_list)<-c("Selected","Feature_plot")
  return(output_list)
}

