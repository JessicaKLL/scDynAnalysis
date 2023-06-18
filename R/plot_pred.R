#' Plot prediction
#' 
#' This function plots the predicted data colouring it different to the original data
#' 
#' @param data The complete data (original + predicted)
#' @param size The width of the plotted line
#' @param base_size Bae font size
#' 
#' @import tidyquant
#' @import ggplot2
#' @import ggformula
#' 
#' @return The plot of the original time steps and the predicted time steps
#' 
#' @export
#' 

plot_pred <- function(data, size = 1, base_size = 14) {
  data$index <- factor(data$index, levels=unique(data$index))
  g <- data %>%
    ggplot(aes(x=index,y=Feature,col=key,group=key)) +
    geom_spline(size=0.75,spar=0.5)+
    theme_tq(base_size = 14) +
    scale_color_tq() + 
    theme(legend.position = "none",axis.text.x = element_blank()) +
    labs(
      title = "",
      x = "Time-points", y = "Expression"
    )
  return(g)
}