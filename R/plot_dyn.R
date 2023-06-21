#' Factor-Feature dynamics
#'
#' This function generates a plot representing the dynamics of the expression of one feature in different patients
#'
#' @param data Data
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

plot_dyn<-function(data){
    p<-ggplot(data,aes(Meta_cells,Gene,group=Individual,col=Individual))+
      geom_point(alpha=0.07)+geom_line(alpha=0.07)+
      geom_spline(size=1,spar = 0.7)
  return(p)
}
