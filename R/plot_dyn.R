#' Factor-Feature dynamics
#'
#' This function generates a plot representing the dynamics of the expression of one feature in different patients
#'
#' @param data List of data splitted by individuals
#' @param cluster Name of the cluster column
#' @param Features The selected features
#'
#' @import ggplot2
#' @import ggformula
#' @import tidyverse
#' @import tibble
#'
#' @return A plot representing the dynamic change of feature and
#' factor expressions.
#'
#' @export
#'

plot_dyn<-function(data,cluster,Features){
  Gene_dyn_Xpatient<-list()
  p<-names(data)
  for (i in 1:length(Features)) {
    if(isTRUE(length(data)==1)){
      ind1<-as.data.frame(data[[i]][,Features[i]])
      ind1$order<-data[[1]][,cluster]
      ind1$individual<-p[1]
      colnames(ind1)<-c("Feature","Meta_cells","Individual")
      df<-ind1
    }
    if(isTRUE(length(data)==2)){
      ind1<-as.data.frame(data[[i]][,Features[i]])
      ind1$order<-data[[2]][,cluster]
      ind1$individual<-p[1]
      colnames(ind1)<-c("Feature","Meta_cells","Individual")
      ind2<-as.data.frame(data[[i]][,Features[i]])
      ind2$order<-data[[2]][,cluster]
      ind2$individual<-p[2]
      colnames(ind2)<-c("Feature","Meta_cells","Individual")
      df<-rbind(ind1,ind2)
    }
    else{
      ind1<-as.data.frame(data[[i]][,Features[i]])
      ind1$order<-data[[1]][,cluster]
      ind1$individual<-p[1]
      colnames(ind1)<-c("Feature","Meta_cells","Individual")
      ind2<-as.data.frame(data[[i]][,Features[i]])
      ind2$order<-data[[2]][,cluster]
      ind2$individual<-p[2]
      colnames(ind2)<-c("Feature","Meta_cells","Individual")
      df<-rbind(ind1,ind2)
      for (i in 3:length(data)) {
        ind<-as.data.frame(data[[i]][,Features[i]])
        ind$order<-data[[i]][,cluster]
        ind$individual<-p[i]
        colnames(ind)<-c("Feature","Meta_cells","Individual")
        df<-rbind(df,ind)
      }
    }
    p<-ggplot(df,aes(Meta_cells,Gene,group=Individual,col=Individual))+
      geom_point(alpha=0.07)+geom_line(alpha=0.07)+
      geom_spline(size=1,spar = 0.7)+
      theme_classic()+theme(axis.text.x = element_blank())+ylab(paste0(Features[i]))+
      labs(title = paste0(Features[i]," dynamics in different individuals"))
    Gene_dyn_Xpatient[[i]]<-p
  }
  names(Gene_dyn_Xpatient)<-Features
  return(Gene_dyn_Xpatient)
}
