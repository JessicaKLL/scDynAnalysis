#' Generate random forest model
#'
#' This function generates a random forest model
#'
#' @param formula Formula describing the model to be fitted
#' @param data The dataframe containing the variables in the model
#' @param split Train Test Split ? (Default=TRUE)
#' @param ntree Number of trees to grow
#' @param strata A variable that is used for stratified sampling (Default = NULL)
#' @param importance Should predictors' importance be assessed?(Default=TRUE,
#' recommended for further analyses)
#' @param proximity Should proimity measure among the rows be calculated? (Default = NULL)
#'
#' @import rsample
#' @import randomForest
#'
#' @return The random forest model with its input data
#'
#' @export


RandomForest<-function(formula,
                       data,
                       split=TRUE,
                       ntree=200,
                       strata=NULL,
                       importance=TRUE,
                       proximity=TRUE){
  if (split) {
    message("Splitting data into training and testing")
    names(data)<-make.names(names(data))
    data_split<-initial_split(data,strata=strata)
    Train_set<-training(data_split)
    Test_set<-testing(data_split)
    model<-list()
    message("Generating the model ...")
    model[[1]]<-randomForest(formula=formula, data = Train_set, ntree=ntree,
                           importance=importance, proximity=proximity)
    model[[2]]<-Train_set
    model[[3]]<-Test_set
    names(model)<-c("RF_model","Train_Set","Test_Set")
    return(model)
  }
  else {
    message("Generating the model ...")
    model<-list()
    names(data)<-make.names(names(data))
    model[[1]]<-randomForest(formula=formula, data = data, ntree=ntree,
                        importance=importance, proximity=proximity)
    model[[2]]<-data
    names(model)<-c("RF_model","Data")
    return(model)
  }
}
