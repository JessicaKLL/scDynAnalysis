#' Factor-Feature dynamics
#'
#' This function generates a plot representing the dynamics of the expression of one feature in different patients
#'
#' @param df Data
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

plot_dyn<-function(df){
    p<-ggplot(df,aes(Meta_cells,Gene,group=Patient,col=Patient))+
      geom_point(alpha=0.07)+geom_line(alpha=0.07)+
      geom_spline(size=1,spar = 0.7)+
      theme_classic()+theme(axis.text.x = element_blank())+ylab(paste0(Features[i]))+
      labs(title = paste0(Features[i]))
    Gene_dyn_Xpatient[[i]]<-p
  }
  return(Gene_dyn_Xpatient)
}
