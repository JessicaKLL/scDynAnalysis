#' Expression heatmap
#'
#' This function plots several heatmaps representing the expression of the genes in
#' different cells at different time-points, taking into account the cell-type.
#'
#' @param data Entire data.frame
#' @param Features The selected features
#'
#' @import pheatmap
#' @import viridis
#'
#'
#' @return A list containing the heatmaps
#'
#' @export
#'

Expr_Heatmap<-function(data, Features){
  data<-split(data, data$cell_type)
  data_au<-rbind(data[[1]],data[[2]])
  for (i in 3:5) {
    data_au<-rbind(data_au,data[[i]])
  }
  data<-data_au
  data<-split(data,data$time_point)
  data_au<-list()
  output_list<-list()
  for (i in 1:length(data)) {
    data_au[[i]]<-data[[i]][,1:length(Features)]
    data_au[[i]]$factor<-data[[i]][,ncol(data[[i]])]
    data_au[[i]]$cell_type<-data[[i]]$cell_type
    n<-length(Features)
    n<-n+1
    x<-t(data_au[[i]][,1:n])
    cluster<-data.frame(rownames(data_au[[i]]),data_au[[i]]$cell_type)
    rownames(cluster)<-cluster[,1]
    cluster[,1]<-NULL
    colnames(cluster)<-"cell_type"
    output_list[[i]]<-pheatmap(x,cluster_rows=F,cluster_cols=F,show_colnames=F,
                               annotation_col=cluster,color=inferno(100),fontsize=12,
                               main=paste0("Feature expression ",names(data)[i]))
  }
  names(output_list)<-names(data)
  return(output_list)
}
