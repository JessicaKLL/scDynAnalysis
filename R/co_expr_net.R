#' Plot co-expression network
#' 
#' This function plots the co-expression network of the features
#' 
#' @param data The feature expression data
#' @param cut_corr Plot features with correlation higher than
#' @param main The title of the plot
#' 
#' @import tidyverse
#' @import tidymodels
#' @import corrr
#' @import ggraph
#' @import igraph
#' 
#' 
#' @return A co-expression network
#' 
#' @export
#' 

co_expr_net<-function(data,cut_corr=0.3,main){
  tidy_cors <- data %>% 
    correlate() %>%
    stretch()
  
  colnames(tidy_cors)<-c("from","to","correlation")
  
  graph_cors <- tidy_cors %>%
    filter(abs(correlation) > cut_correlation) %>%
    graph_from_data_frame(directed = FALSE)
  
  cor_net<-ggraph(graph_cors) +
    geom_edge_link(aes(edge_alpha = abs(correlation), color = correlation)) +
    guides(edge_alpha = "none", edge_width = "none") +
    scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) +
    geom_node_point(color = "white", size = 5) +
    geom_node_text(aes(label = name), repel = TRUE)+
    labs(title = main)

  return(cor_net)
}
