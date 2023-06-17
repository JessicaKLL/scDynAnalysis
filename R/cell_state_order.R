#' Ordering the cell states
#' 
#' This function orders the cell states inside one time point of the biological process based on the 
#' distances between cell states
#' 
#' @param data Your data (one time point)
#' @param root The starting cell state of the ordering process
#' @param cluster Name of the cluster column
#' @param feat_num Number of features
#' @param dist_method The distance method used to generate the distance matrix
#' 
#' 
#' @import dplyr
#' @import cluster
#' @import tidyverse
#' 
#' @return A 
#' 

cell_state_order<-function(data,root,cluster,feat_num,dist_method="minkowski"){
  cell_state<-c(root)
  d<-dist(data[,1:feat_num],method = dist_method)
  d<-as.matrix(d)
  rownames(d)<-data[,cluster]
  colnames(d)<-data[,cluster]
  d<-as.dist(d)
  d<-as.matrix(d)
  d_au<-d[,root]
  d<-d[-which(rownames(d) %in% cell_state),]
  d_au<-as.data.frame(d_au)
  d_au$clust<-rownames(d_au)
  d_au<-as.data.frame(d_au[order(d_au$d_au),])
  new_root<-d_au$clust[2]
  cell_state<-append(cell_state,new_root)
  i<-2
  while(i!=nrow(data)){
    d_au<-d[,new_root]
    d<-d[-which(rownames(d) %in% cell_state),]
    d_au<-as.data.frame(d_au)
    d_au$clust<-rownames(d_au)
    d_au<-as.data.frame(d_au[order(d_au$d_au),])
    new_root<-d_au$clust[2]
    cell_state<-append(cell_state,new_root)
    i<-i+1
  }
  return(cell_state)
}