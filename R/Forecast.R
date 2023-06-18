#' Forecasting future time steps
#' 
#' This function forecasts the future time steps of a selected feature
#' 
#' @param data Your time step dataframe
#' @param batch_size The number of training examples in one forward/backward pass of a 
#' RNN before a weight update. It must be evenly divisible into both the training and testing lengths
#' @param train_length Number of training time steps
#' @param tsteps The number of lags included in the training/testing set
#' @param epochs The total number of forward/backward pass iterations
#' 
#' @import glue
#' @import forcats
#' @import timetk
#' @import tidyquant
#' @import tibbletime
#' @import cowplot
#' @import recipes
#' @import rsample
#' @import yardstick
#' @import keras
#' 
#' @return A dataframe with both input and forecasted time steps
#' 
#' @export
#' 

Forecast <- function(data,batch_size=25,train_length=150,
                                      tsteps=1,epochs=500) {
  df <- data
  
  colnames(df)<-c("index","Feature")
  
  lag_train_tbl <- df %>%
    mutate(Feature_lag = lag(Feature, n = lag_setting)) %>%
    filter(!is.na(Feature_lag)) %>%
    tail(train_length)
  
  x_train_vec <- lag_train_tbl$Feature_lag
  x_train_arr <- array(data = x_train_vec, dim = c(length(x_train_vec), 1, 1))
  
  y_train_vec <- lag_train_tbl$Feature
  y_train_arr <- array(data = y_train_vec, dim = c(length(y_train_vec), 1))
  
  x_test_vec <- y_train_vec %>% tail(lag_setting)
  x_test_arr <- array(data = x_test_vec, dim = c(length(x_test_vec), 1, 1))
  
  # LSTM Model
  model <- keras_model_sequential()
  model %>%
    layer_lstm(units            = 20, 
               input_shape      = c(tsteps, 1), 
               batch_size       = batch_size,
               return_sequences = TRUE, 
               stateful         = TRUE,
               activation = "relu") %>% 
    layer_dropout(rate = 0.2) %>% 
    layer_lstm(units            = 50, 
               return_sequences = TRUE, 
               stateful         = TRUE) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_lstm(units            = 50, 
               return_sequences = TRUE, 
               stateful         = TRUE) %>% 
    layer_dropout(rate = 0.1) %>%
    layer_lstm(units            = 20, 
               return_sequences = FALSE, 
               stateful         = TRUE) %>% 
    layer_dense(units = 5)
  
  model %>% 
    compile(loss = 'mae', optimizer = 'adam')
  
  # Fitting LSTM
  for (i in 1:epochs) {
    model %>% fit(x          = x_train_arr, 
                  y          = y_train_arr, 
                  batch_size = batch_size,
                  epochs     = 1, 
                  verbose    = 1, 
                  shuffle    = FALSE)
    
    model %>% reset_states()
    cat("Epoch: ", i)
    
  }
  
  # Predict and Return Tidy Data (MODIFIED)
  # Make Predictions
  pred_out <- model %>% 
    predict(x_test_arr, batch_size = batch_size) %>%
    .[,1]
  
  # Make future index 
  idx<-c()
  n<-length(pred_out)
  for (i in 1:n) {
    x<-paste0("PredTP_",i)
    idx<-append(idx,x)
  }
  
  # Transform Features
  df <- tibble(
    index=df$index,
    Feature=df$Feature
  )
  
  pred_tbl <- tibble(
    index   = lag_test_tbl$index,
    Feature   = pred_out
  )
  
  # Combine actual data with predictions
  tbl_1 <- df %>%
    add_column(key = "actual")
  tbl_3 <- pred_tbl %>%
    add_column(key = "predict")
  
  colnames(tbl_1)<-colnames(tbl_3)
  
  complete_tbl<-rbind(tbl_1,tbl_3)
  
  return(complete_tbl)
}