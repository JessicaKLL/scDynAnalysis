#' Plot extracted decision tree from forest
#'
#' This function plots a decision tree extracted from a random forest model.
#'
#'
#' @param model The random forest model
#' @param tree_num The index of the tree which you want to extract from the model
#' @param main The title of the plot
#'
#' @import randomForest
#' @import dplyr
#' @import ggraph
#' @import igraph
#'
#' @return A decision tree
#'
#' @export
#'

rf_tree_func <- function(model,
                      tree_num, main) {

  # get tree by index
  tree <- randomForest::getTree(model,
                                k = tree_num,
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))

  # prepare data frame for graph
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))
  # convert to graph and delete the last node that we don't want to plot
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")
  # set node labels
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$node_label
  V(graph)$leaf_label <- as.character(tree$prediction)
  V(graph)$split <- as.character(round(tree$`split point`, digits = 3))

  # plot
  plot <- ggraph(graph, "tree", circular=F) +
    theme_bw() +
    geom_edge_elbow() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE, size=3, fontface = "bold") +
    #geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label, alpha=0.1), na.rm = TRUE,
                    repel = TRUE, colour = "black", fontface = "bold", show.legend = FALSE, size=2.5) +
    labs(title = main)+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 15))
  return(plot)
}
