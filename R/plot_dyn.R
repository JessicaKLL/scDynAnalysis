#' Factor-Feature dynamics
#'
#' This function generates a plot representing the dynamics of the expression of one feature in different patients
#'
#' @param data Data
#' @param spar Smoothing parameter
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

plot_dyn<-function(data,spar=0.5){
    p<-ggplot(data,aes(Meta_cells,Gene,group=Individual,col=Individual))+
      geom_point(alpha=0.07)+geom_line(alpha=0.07)+
      geom_spline(size=1,spar = spar)
  return(p)
}
