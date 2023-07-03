#' Expression heatmap
#'
#' This function plots a heatmaps representing the expression of the genes in
#' different cells at different time-points.
#'
#' @param data Entire data.frame
#' @param Features The selected features
#' @param main Title of the plot
#'
#' @import pheatmap
#' @import viridis
#'
#'
#' @return A heatmap
#'
#' @export
#'

Expr_Heatmap<-function(data, Features,main=""){
    data_au<-data[,which(colnames(data) %in% Features)]
    data_au$time_point<-data$time_point
    n<-length(Features)
    x<-t(data_au[,1:n])
    cluster<-data.frame(rownames(data_au),data_au$time_point)
    rownames(cluster)<-cluster[,1]
    cluster[,1]<-NULL
    colnames(cluster)<-"time_point"
    breaksList = seq(0, 100, by = 1)
    output<-pheatmap(x,cluster_rows=F,cluster_cols=F,show_colnames=F,show_rownames = T,
                     annotation_col=cluster,color=turbo(length(breaksList)),
                     breaks=breaksList,fontsize=12,main = main)
    return(output)
}
