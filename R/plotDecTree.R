#' Plot the decision tree of the model
#'
#' This function plots a decision tree, it can be extracted from the random forest model, or
#' from generating a decision tree model.
#'
#' @param model_type Select which model do you want to use (rf=randomForest, tree=tree model)
#' @param generate Should a
#' @param model The model. If no model provided, a selected (rf or tree) will be generated
#' @param tree_num The index of the tree which you want to extract from the model. Only when 'model_type="rf"'
#' @param num_tree The number of trees you want to generate.
#' Only when 'model_type="rf"' and no model is provided (Default=200)
#' @param main The title of the plot
#' @param data The data that will be used to generate the model. Only when no model is provided
#' @param formula The formula to generate the model. Only when no model is provided
#' @param type The type of plot that you want to generate from the tree model. Only when 'model_type="tree"'
#'
#' @import rpart
#' @import rpart.plot
#' @import randomForest
#' @import dplyr
#' @import ggraph
#' @import igraph
#'
#' @return A decision tree
#'
#' @export
#'

plotDecTree<-function(model_type=c("rf","tree"),model,tree_num=2,num_tree=200,main="Decision tree",data,formula,type=0){
  if (model_type=="rf"){
    if (missing(model)) {
      model<-RandomForest(formula,data,split=FALSE,ntree=num_tree,strata=NULL,importance=T,proximity=F)
      plot<-rf_tree_func(model[[1]],tree_num, main)
      return(plot)
    }
    else{
      plot<-rf_tree_func(model,tree_num, main)
      return(plot)
    }
  }
  if (model_type=="tree"){
    if (missing(model)) {
      x<-rpart(formula ,data)
      rpart.plot(x,type = 0,main=main)
      plot<-recordPlot()
      return(plot)
    }
    else{
      rpart.plot(model,type = 0,main=main)
      plot<-recordPlot()
      return(plot)
    }
  }
}
