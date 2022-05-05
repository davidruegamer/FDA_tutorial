library(tensorflow)

# function to use for applying a pre-trained net
pretrained_fitting <- function(
  training_data_x,
  training_data_y,
  validation_data,
  width = 150, height = 150,
  pretrained_model = "imagenet",
  top_net = 
    function(x){
      x %>% layer_flatten() %>%
        layer_dense (units = 256, kernel_regularizer = regularizer_l2(0.01)) %>%
        layer_activation(activation = "relu") %>%
        layer_batch_normalization() %>%
        layer_dropout(rate = 0.5) %>%
        layer_dense (units = 101, kernel_regularizer = regularizer_l2(0.01)) %>%
        layer_activation(activation = "relu") %>%
        layer_batch_normalization() %>%
        layer_dropout(rate = 0.25) %>%
        layer_dense (units = 101, activation = "linear")
    },
  epochs = 200,
  batch_size = 16,
  callbacks = list(
    callback_early_stopping(patience = 20, restore_best_weights = TRUE),
    callback_learning_rate_scheduler(function(epoch, lr){
      if(epoch<10){
        return(lr)
      }else{
        return(lr * tf$math$exp(-0.1))
      }
    })
    # callback_learning_rate_scheduler(
    #   tf$keras$experimental$CosineDecayRestarts(.02, 10, t_mul = 2, m_mul = .7)
    # )
  ),
  optimizer = "rmsprop",
  loss = "mse",
  verbose = TRUE
)
{
  
  library(keras)
  
  conv_base <- application_densenet121(
    weights = pretrained_model,
    include_top = FALSE,
    input_shape = c(width, height, 3)
  )
  
  # Create model ------------------------------------------------------------
  model <- keras_model_sequential() %>%
    conv_base %>%
    top_net()
  
  # Consider freezing bottom weights to prevent overfit
  freeze_weights(conv_base)
  
  # run model ------------------------------------------------------------
  model %>% compile(
    optimizer = optimizer,
    loss = loss
  )
  
  history <- model %>% fit(
    verbose = verbose,
    view_metrics = FALSE,
    training_data_x,
    training_data_y,
    epochs = epochs,
    batch_size = batch_size,
    validation_data = validation_data,
    callbacks = callbacks
  )
  
  return(model)
  
}

# CNN Block
cnn_block <- function(filters, kernel_size, pool_size, rate, input_shape = NULL){
  function(x){
    x %>%
      layer_conv_2d(filters, kernel_size, padding="same", input_shape = input_shape) %>%
      layer_activation(activation = "relu") %>%
      layer_batch_normalization() %>%
      layer_max_pooling_2d(pool_size = pool_size) %>%
      layer_dropout(rate = rate)
  }
}

# Conv NN
conv_fitting <- function(
  training_data_x,
  training_data_y,
  validation_data,
  width = 150, height = 150,
  filters = c(32, 64, 128),
  kernel_sizes = c(3, 3, 3),
  pool_sizes = c(3, 2, 2),
  dropout_rates = c(0.25, 0.25, 0.25),
  top_layer = function(x){
    x %>% 
      layer_flatten() %>%
      layer_dense(256) %>%
      layer_activation(activation = "relu") %>%
      layer_batch_normalization() %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(101)
  },
  epochs = 200,
  batch_size = 16,
  callbacks = list(
    # callback_early_stopping(patience = 25, restore_best_weights = TRUE),
    callback_learning_rate_scheduler(
      tf$keras$experimental$CosineDecayRestarts(.02, 10, t_mul = 2, m_mul = .7)
    )
  ),
  optimizer = "rmsprop",
  loss = "mse",
  verbose = TRUE
){
  
  cnn_blocks <- lapply(1:length(filters), function(i){
    
    inp_sh <- NULL
    if(i==1) inp_sh <- shape(width, height, 3)
    cnn_block(filters = filters[i], 
              kernel_size = rep(kernel_sizes[i],2), 
              pool_size = rep(pool_sizes[i],2), 
              rate = dropout_rates[i],
              inp_sh)
    
  })

    # Create model --------------------------------------------------------------
    
    model <- keras_model_sequential()
    
    for(i in 1:length(cnn_blocks)){
      
      model <- model %>% cnn_blocks[[i]]()
      
    }
    
    model <- model %>% top_layer()

    model %>% compile(
      optimizer = optimizer,
      loss = loss
    )
    
    
    history <- model %>% fit(
      verbose = verbose,
      view_metrics = FALSE,
      training_data_x,
      training_data_y,
      epochs = epochs,
      batch_size = batch_size,
      validation_data = validation_data,
      callbacks = callbacks
    )
    
    return(model)
    
}

funvar_to_image <- function(data,
                            wid = 150,
                            ht = 150)
{
  
  library (tidyverse)
  library (imager)
  
  train_mat <- data[map_lgl (data, is.matrix)]
  train_x <- train_mat[grepl ("accl|vel|angle", names (train_mat))]
  axes <- c("ap", "ml", "vt")
  prednames <- unique (str_remove_all(pred_vars, "_ap|_ml|_vt"))
  outnames <- unique (str_remove_all(response_vars, "_ap|_ml|_vt"))
  
  ### Restructure predictors
  
  temp_array <- array (dim = c(nrow (train_x[[1]]), ncol (train_x[[1]]), 9))
  train_x_array <- array (dim = c(nrow (train_x[[1]]), ncol (train_x[[1]]), 9, 3))
  
  for (n in seq_along (axes)) {
    
    temp <- train_x[grepl (axes[n], names (train_x))]
    
    for (m in seq_along (temp)) {
      
      temp_array[, , m] <- temp[[m]]
      
    }
    
    train_x_array[, , , n] <- temp_array
    
  }
  
  resize(train_x_array, size_y = wid, size_z = ht, interpolation_type = 5) %>% as.array()

}
