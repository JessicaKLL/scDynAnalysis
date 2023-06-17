#' Divisive hierarchical clustering
#' 
#' This function applies a hierarchical clustering to each time point of your data,
#' generating cluster considered as cell states.
#' 
#' @param data Your data (must be splitted by time points)
#' @param k Number of clusters (default=200)
#' @param dist_method The method to generate the distance matrix for clustering (default="minkowski")
#' @param clust_method The method to clusterize the data (default="ward.D2")
#' 
#' @return The original data with a column of clusters
#' 
#' @import cluster

divHier<-function(data,k=200,dist_method="minkowski",clust_method="ward.D2"){
  tp_hclust<-list()
  tp<-names(data)
  for (x in 1:length(data)) {
    d<-dist(data[[x]],method = dist_method)
    final_clust<-hclust(d,method = clust_method)
    groups<-cutree(final_clust,k=k)
    tp_hclust[[x]]<-data.frame(groups)
  }
  
  names(tp_hclust)<-tp
  
  for (x in 1:length(tp_hclust)) {
    for (i in 1:nrow(tp_hclust[[x]])) {
      tp_hclust[[x]]$groups[i]<-paste0(tp[x],"_",tp_hclust[[x]]$groups[i])
    }
  }
  
  for (x in 1:length(data)) {
    data[[x]]$clusters<-tp_hclust[[x]]$groups
  }
  
  return(data)
}
