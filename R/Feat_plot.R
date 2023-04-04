#' Feature plot
#'
#' This function generates a list containing the feature plot of each selected feature.
#' Feature plot: Coloring single cells on a dimensional reduction plot according to a 'feature'
#'
#' @param object The Seurat object
#' @param Features Vector containing the names of the selected features
#' @param raster Convert points to raster format (Default=NULL)
#' @param cells Select cells to be plotted (Default=NULL)
#' @param order Boolean determining whether to plot cells in order of expression (Default=NULL)
#'
#' @import Seurat
#'
#' @return List containing the feature plot of each selected feature
#'
#' @export

feat_plot<-function(object,Features,raster=NULL,cells=NULL,order=NULL){
  out_list<-list()
  for (i in 1:length(Features)) {
    out_list[[i]]<-FeaturePlot(CITE_data, features = Features[i],raster=raster,
                                     cells = cells,order=order)
  }
  return(out_list)
}
