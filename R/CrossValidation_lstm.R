#' Cross-Validation of lstm model
#' 
#' This function runs a cross-validation lstm model to check if your parameters are 
#' appropriate for the forecasting of future trends of the feature
#' 
#' @param data Your time step data
#' @param Feature Interested feature
#' @param periods_train Number of training time steps
#' @param periods_test Number of testing time steps
#' @param skip_span Distribute samples into n sets that span the entire historical data
#' @param batch_size The number of training examples in one forward/backward pass of a 
#' RNN before a weight update. It must be evenly divisible into both the training and testing lengths
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
#' @return A dataframe with both original and predicted testing time steps
#' 
#' @export
#' 

CrossValidation_lstm<-function(data,Feature,periods_train=300,periods_test=100,
                                skip_span=10,batch_size=50,tsteps=1,epochs=1000){
  
  tk<-Generate_time_series_featDF(data,Feature = Feature)
  
  rolling_origin_resamples <- rolling_origin(
    tk,
    initial    = periods_train,
    assess     = periods_test,
    cumulative = ,
    skip       = skip_span
  )
  split<-rolling_origin_resamples$splits[[1]]
  df_trn <- training(split)
  colnames(df_trn)<-c("index","Feature")
  df_trn_au<-as.data.frame(df_trn$index)
  df_trn<-cbind(df_trn_au,df_trn$Feature)
  colnames(df_trn)<-c("index","Feature")
  df_tst <- testing(split)
  colnames(df_tst)<-c("index","Feature")
  df_tst_au<-as.data.frame(df_tst$index)
  df_tst<-cbind(df_tst_au,df_tst$Feature)
  colnames(df_tst)<-c("index","Feature")
  df <- bind_rows(
    df_trn %>% add_column(key = "training"),
    df_tst %>% add_column(key = "testing")
  ) 
  
  lag_setting  <- nrow(df_tst)
  train_length <- nrow(df_trn)
  
  # Training Set
  lag_train_tbl <- df %>%
    mutate(Feature_lag = lag(Feature, n = lag_setting)) %>%
    filter(!is.na(Feature_lag)) %>%
    filter(key == "training") %>%
    tail(train_length)
  
  x_train_vec <- lag_train_tbl$Feature_lag
  x_train_arr <- array(data = x_train_vec, dim = c(length(x_train_vec), 1, 1))
  #x_train_arr <- array(data = x_train_vec, dim = c(train_length, lag_setting, 1))
  
  
  y_train_vec <- lag_train_tbl$Feature
  y_train_arr <- array(data = y_train_vec, dim = c(length(y_train_vec), 1))
  
  # Testing Set
  lag_test_tbl <- df%>%
    mutate(
      Feature_lag = lag(Feature, n = lag_setting)
    ) %>%
    filter(!is.na(Feature_lag)) %>%
    filter(key == "testing")
  
  x_test_vec <- lag_test_tbl$Feature_lag
  x_test_arr <- array(data = x_test_vec, dim = c(length(x_test_vec), 1, 1))
  #x_test_arr <- array(data = x_test_vec, dim = c(nrow(lag_test_tbl), lag_setting, 1))
  
  y_test_vec <- lag_test_tbl$Feature
  y_test_arr <- array(data = y_test_vec, dim = c(length(y_test_vec), 1))
  
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
    compile(loss= 'mae', optimizer_adam(learning_rate=0.001))
  model
  
  
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
  
  
  # Make Predictions
  pred_out <- model %>% 
    predict(x_test_arr, batch_size = batch_size) 
  
  # Retransform values
  
  pred_tbl <- tibble(
    index   = lag_test_tbl$index,
    Feature   = pred_out
  )
  
  df_trn<-tibble(
    index = df_trn$index,
    Feature = df_trn$Feature
  )
  
  df_tst<-tibble(
    index=df_tst$index,
    Feature = df_tst$Feature
  )

  tbl_1 <- df_trn %>%
    add_column(key = "actual")
  
  tbl_2 <- df_tst %>%
    add_column(key = "actual")
  
  tbl_3 <- pred_tbl %>%
    add_column(key = "predict")
  
  colnames(tbl_1)<-colnames(tbl_3)
  colnames(tbl_2)<-colnames(tbl_3)
  complete_tbl<-rbind(tbl_1,tbl_2)
  complete_tbl<-rbind(complete_tbl, tbl_3)
  return(complete_tbl)
}
