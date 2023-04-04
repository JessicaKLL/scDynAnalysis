#' Check quality
#'
#' This function checks the quality of the model
#'
#' @param model The random forest model
#' @param test_set Data which variable/factor/output is predicted by the model
#' @param real_output The real variable/factor/output of the test data
#'
#' @import caret
#' @import dplyr
#'
#' @return List containing the predicted data, the confusion matrix, the model quality table and a heatmap representing
#' the correctness of the prediction.
#'
#' @export

model_quality<-function(model,test_set,real_output){
  predicted_set<-predict(model,newdata=test_set,type="response")
  predicted_set<-factor(predicted_set,levels(real_output))
  real<-as.data.frame(table(real_output))
  m<-confusionMatrix(predicted_set, real_output)
  confusion_matrix<-m$table
  Quality<-m$byClass
  confusion<-data.frame(m$table)
  confusion$perc<-confusion$Freq/real$Freq*100
  colnames(confusion)<-c("Predicted","Real","Count","Perc")
  plot<- confusion %>%
    ggplot(aes(x=Real,y=Predicted,fill=Perc)) +
    geom_tile()+scale_fill_viridis_c() + theme(axis.text.x = element_text(angle = 90)) +
    geom_text(aes(label=round(Perc,3)))+
    labs(title = "Heatmap comparing predicted and real values")
  output<-list()
  output[[1]]<-predicted_set
  output[[2]]<-confusion_matrix
  output[[3]]<-Quality
  output[[4]]<-plot
  names(output)<-c("Prediction","Confusion_Matrix","Model_Quality","Heatmap")
  return(output)
}
