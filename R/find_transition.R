#' Find transition points between two time points
#' 
#' This function finds out the transition points between two time points based on their correlation
#' 
#' @param data You data
#' @param cluster Name of the cluster column
#' @param feat_num Number of features
#' 
#' @import corrr
#' 
#' @return A list that contains a correlation matrix and the two clusters with highest correlation, 
#' considered as transition points between the two time points.
#' 

find_transition<-function(data,cluster,feat_num){
  n<-split(data,data$time_point)
  n<-rbind(n$IP,n$Peak)
  CL_G<-data.frame(n[,1:feat_num])
  CL_G<-as.data.frame(t(CL_G))
  colnames(CL_G)<-n[,cluster]
  corr_cl<-correlate(CL_G)
  corr_cl<-stretch(corr_cl)
  
  stop<-nrow(n)*nrow(n)/2
  while(stop!=nrow(corr_cl)){
    for (i in 1:nrow(corr_cl)) {
      if(isTRUE(strsplit(corr_cl$x[i], "[_]")[[1]][1] == strsplit(corr_cl$y[i], "[_]")[[1]][1])){
        corr_cl<-corr_cl[-i,]
      }
    }
  }
  
  identical <- corr_cl[which.max(corr_cl$r),]
  
  output<-list()
  
  output[[1]]<-corr_cl
  output[[2]]<-identical
  names(output)<-c("Corr_Mat","Identical")
  
  return(output)
}
