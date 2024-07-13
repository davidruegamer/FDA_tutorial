### Section 6: Machine Learning Approaches ###

# The following code is only reproducible using
# TensorFlow >= 2.10 and Keras >= 2.10, but Keras < 3.
# Due to the lack in support of earlier versions of 
# these packages on Windows, the code might further
# only be executable on Linux-based systems.

# library import
library(refund) # functional regression models for comparison
library(FuncNN) # neural networks with functional input
library(FDboost) # Boosting functional regression
library(tidyverse) # data wrangling
library(ggplot2) # plotting
source("Rcode/nn_helpers.R") # neural networks

# create folder
if(!dir.exists("results"))
  dir.create("results")

# data import
dta <- readRDS("data/data_comb.RDS")
names(dta)

### 6.4 Comparison ###

# Get variables---------------------------------------------------------------

response_vars <- names(dta)[grep("grf|knee_moment|hip_moment|ankle_moment", 
                                 names(dta))]
pred_vars <- names(dta)[grep("angle|vel|accl", names(dta))]

# Restructure data -----------------------------------------------------------

dta_mat <- dta[map_lgl (dta, is.matrix)]
dta_x <- dta_mat[grepl ("accl|vel|angle", names (dta_mat))]
dta_y <- dta_mat[grepl ("grf|knee_moment|hip_moment|ankle_moment", 
                        names(dta_mat))]
axes <- c("ap", "ml", "vt")
prednames <- unique (str_remove_all(pred_vars, "_ap|_ml|_vt"))
outnames <- unique (str_remove_all(response_vars, "_ap|_ml|_vt"))

set.seed(42)

train_ind <- sample(1:nrow(dta_y[[1]]), floor(0.8*nrow(dta_y[[1]])))
test_ind <- setdiff(1:nrow(dta_y[[1]]), train_ind)

# define final data for comparison
x_sca <- model.matrix(~ -1 + ., data = as.data.frame(
  dta[c("study", "sex", "cond", "side", "step", "age", "ht", "wt")]
))

y_train <- dta_y$ankle_moment_ap[train_ind,]
x_train_fun <- lapply(dta_x, function(x) x[train_ind,])
x_train_sca <- x_sca[train_ind,]

y_test <- dta_y$ankle_moment_ap[test_ind,]
x_test_fun <- lapply(dta_x, function(x) x[test_ind,])
x_test_sca <- x_sca[test_ind,]

# FDA as image
data_array <- funvar_to_image(dta)
x_train_fun_array <- data_array[train_ind,,,]
x_test_fun_array <- data_array[test_ind,,,]

# Another split for validation
train_train_ind <- sample(1:nrow(x_train_fun[[1]]), 
                          floor(0.8*nrow(x_train_fun[[1]])))
train_val_ind <- setdiff(1:nrow(x_train_fun[[1]]), train_train_ind)

# for pffr
cycle <- dta$cycle
train <- dta[names(dta)!="cycle"]
train <- lapply(train, function(x) if(is.null(dim(x))) x[train_ind] else 
  x[train_ind,])
test <- dta[names(dta)!="cycle"]
test <- lapply(test, function(x) if(is.null(dim(x))) x[test_ind] else 
  x[test_ind,])

# Models ---------------------------------------------------------------------

################ FuncNN #################
fit_funcNN = fnn.fit(resp = y_train,
                     func_cov = x_train_fun,
                     scalar_cov = x_train_sca,
                     hidden_layers = 6,
                     neurons_per_layer = c(24, 24, 24, 24, 24, 58),
                     activations_in_layers = c("relu", "relu", "relu", 
                                               "relu", "relu", "linear"),
                     domain_range = list(c(1, 101)),
                     learn_rate = 0.001,
                     epochs = 100,
                     raw_data = TRUE,
                     early_stopping = TRUE)
# 
# # Running prediction, gets probabilities
prediction_funcNN = fnn.predict(fit_funcNN,
                                func_cov = x_test_fun,
                                scalar_cov = x_test_sca,
                                domain_range = list(c(1, 101)),
                                raw_data = TRUE)

saveRDS(prediction_funcNN, file="results/prediction_FuncNN.RDS")

rm(fit_funcNN, prediction_funcNN); gc()

################ Pre-trained Deep Network ##################

fit_imagenet <- pretrained_fitting(
  x_train_fun_array[train_train_ind,,,],
  y_train[train_train_ind,],
  list(x_train_fun_array[train_val_ind,,,], 
       y_train[train_val_ind,])
)

prediction_imagenet <- fit_imagenet %>% predict(x_test_fun_array)

saveRDS(prediction_imagenet, file="results/prediction_ImageNet.RDS")

rm(fit_imagenet, prediction_imagenet); gc()

################ Convolutional Neural Network ##################

fit_cnn <- conv_fitting(
  x_train_fun_array[train_train_ind,,,],
  y_train[train_train_ind,],
  list(x_train_fun_array[train_val_ind,,,], 
       y_train[train_val_ind,])
)

prediction_cnn <- fit_cnn %>% predict(x_test_fun_array)

saveRDS(prediction_cnn, file="results/prediction_CNN.RDS")

rm(fit_cnn, prediction_cnn); gc()

################ PFFR ####################

response <- response_vars[1]

# Create formula
form <- paste(response, " ~ 1 + ", paste(
  paste0("ff(", pred_vars,
         ", yind=cycle, xind=cycle)"),
  collapse = " + "),
  " + age",
  " + ht",
  " + wt",
  " + sex",
  "+ s(cond, bs = 're')",
  "+ s(id, bs = 're')"
)

train <- as.data.frame(lapply(train, function(x) if(!is.null(dim(x))) 
  return(I(x)) else x))
train$id <- as.factor(train$id)
train$cond <- as.factor(train$cond)
train$sex <- as.factor(train$sex)
test <- as.data.frame(lapply(test, function(x) if(!is.null(dim(x))) 
  return(I(x)) else x))
test$id <- as.factor(test$id)
test$cond <- as.factor(test$cond)
test$sex <- as.factor(test$sex)

# initialize the model
m <- pffr(as.formula(form),
          yind = 1:101,
          algorithm = "bam",
          # discrete = TRUE,
          data = train)

prediction_pffr <- m %>% predict(test)

saveRDS(prediction_pffr, file="results/prediction_pffr.RDS")

rm(m, prediction_pffr); gc()

################ FDboost ####################

response <- response_vars[1]

# Create formula
form <- paste(response, " ~ 1 + ", paste(
  paste0("bsignal(", pred_vars,
         ", cycle)"),
  collapse = " + "),
  " + bbsc(age)",
  " + bbsc(ht)",
  " + bbsc(wt)",
  " + bolsc(sex, df = 2)",
  "+ brandomc(cond)"
)

train <- as.list(train)
train$cycle <- cycle
train[pred_vars] <- lapply(train[pred_vars], function(x) scale(x, scale=F))
test <- as.list(test)
test$cycle <- cycle
test[pred_vars] <- lapply(test[pred_vars], function(x) scale(x, scale=F))

# initialize the model
m <- FDboost(as.formula(form),
             data = train,
             timeformula = ~ bbs(cycle, df = 5),
             control=boost_control(mstop = 1000, nu = 0.1))

set.seed(123)
appl1 <- applyFolds(m, folds = mboost::cv(rep(1, length(unique(m$id))), 
                                          B = 5), 
                    grid = 1:1000)
## plot(appl1)
m[mstop(m)]

prediction_fdboost <- m %>% predict(test)

saveRDS(prediction_fdboost, file="results/prediction_FDboost.RDS")

rm(m, prediction_fdboost); gc()

################ Functional Intercept ####################

response <- response_vars[1]

# Create formula
form <- paste(response, " ~ 1")

# initialize the model
mint <- pffr(as.formula(form),
             yind = 1:101,
             algorithm = "bam",
             data = train)

prediction_intercept <- mint %>% predict(test[c("id")])

saveRDS(prediction_intercept, file="results/prediction_intercept.RDS")

rm(mint, prediction_intercept); gc()

###############################################################
######################### Comparison ##########################
###############################################################

results <- c(list(y_test), lapply(list.files("results", full.names = T), 
                                  function(x) as.matrix(readRDS(x))))
nams <- c("truth", gsub("prediction_(.*)\\.RDS", "\\1", list.files("results")))
resultsDF <- do.call("rbind", lapply(1:length(results), function(i) 
  data.frame(value = c(results[[i]]), time = rep(1:101, each = nrow(y_test)),
             obs = rep(1:nrow(y_test), 101), what = nams[i])))

random_g <- mean(sapply(1:100, function(i) 
  mean(sqrt(apply((y_train[sample(1:nrow(y_train), nrow(y_test)),]-
                     results[[1]])^2, 1, sum)))))

rmseDF <- data.frame(rmse = paste0("RMSE: ", round(
  sapply(2:length(results), function(i) 
    mean(sqrt(apply((results[[i]]-results[[1]])^2, 1, sum)))), 4)),
  what = nams[-1], obs = 1)

meth_lab <- c("CNN", "FDboost", "FuncNN", "ImageNet", "Intercept", "pffr")

rmseDF$what <- factor(rmseDF$what, levels = unique(rmseDF$what),
                      labels = meth_lab)
resultsDF$what <- factor(resultsDF$what, levels = unique(resultsDF$what),
                         labels = c("Truth", meth_lab))

# Filter the results
filtered_results <- subset(resultsDF, !what %in% c("Intercept"))
filtered_rmse <- subset(rmseDF, !what %in% c("Intercept"))

# Plotting each method in a separate plot
palette <- RColorBrewer::brewer.pal(8, "Dark2")
palette <- col2rgb(palette)
palette <- rgb(palette[1,] / 255, palette[2,] / 255, palette[3,] / 255, alpha = 0.25)

# Figure 13

par(mfrow = c(2,3), cex=0.8)
for (method in unique(filtered_results$what)) {
  single_result <- subset(filtered_results, what == method)
  matplot(t(matrix(single_result$value, ncol=101)), 
          x = seq(0,1,l=101), 
          col = palette[as.integer(single_result$what)], 
          type = 'l', ylim = c(-1.1,1.1), 
          xlab = 'relative time', ylab = 'value', main = method, bty = "n")
  
  single_rmse <- subset(filtered_rmse, what == method)
  if(method!="Truth") text(x = min(single_result$time)/101, y = -0.82, 
                           labels = single_rmse$rmse, pos = 4)
}

