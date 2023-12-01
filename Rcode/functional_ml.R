setwd("..")
rm(list=ls())
gc()
# Load packages --------------------------------------------------------------
library(refund) # functional regression models for comparison
library(FuncNN) # neural networks with functional input
# load newest version from Github:
# https://github.com/gpfda/GPFDA-dev/issues/2#issuecomment-1484055581
library(GPFDA) # GPs for FDA
library(FDboost) # Boosting functional regression
library(FunFor) # functional random forest
library(xgboost) # for implicit ML approach
source("Rcode/nn_helpers.R") # neural networks

# data wrangling
library(tidyverse)
# plotting
library(ggplot2)

# data import-----------------------------------------------------------------

dta <- readRDS("data/data_comb.RDS")

# Get variables---------------------------------------------------------------

response_vars <- names(dta)[grep("grf|knee_moment|hip_moment|ankle_moment", names(dta))]
pred_vars <- names(dta)[grep("angle|vel|accl", names(dta))]

# Restructure data -----------------------------------------------------------

dta_mat <- dta[map_lgl (dta, is.matrix)]
dta_x <- dta_mat[grepl ("accl|vel|angle", names (dta_mat))]
dta_y <- dta_mat[grepl ("grf|knee_moment|hip_moment|ankle_moment", names (dta_mat))]
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
train_train_ind <- sample(1:nrow(x_train_fun[[1]]), floor(0.8*nrow(x_train_fun[[1]])))
train_val_ind <- setdiff(1:nrow(x_train_fun[[1]]), train_train_ind)

# for pffr
cycle <- dta$cycle
train <- dta[names(dta)!="cycle"]
train <- lapply(train, function(x) if(is.null(dim(x))) x[train_ind] else x[train_ind,])
test <- dta[names(dta)!="cycle"]
test <- lapply(test, function(x) if(is.null(dim(x))) x[test_ind] else x[test_ind,])

# Models ---------------------------------------------------------------------

################ FuncNN #################
fit_funcNN = fnn.fit(resp = y_train,
                    func_cov = x_train_fun,
                    scalar_cov = x_train_sca,
                    hidden_layers = 6,
                    neurons_per_layer = c(24, 24, 24, 24, 24, 58),
                    activations_in_layers = c("relu", "relu", "relu", "relu", "relu", "linear"),
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

saveRDS(prediction_funcNN, file="results/prediction_funcNN.RDS")

rm(fit_funcNN, prediction_funcNN); gc()

################ Pre-trained Deep Network ##################

fit_imagenet <- pretrained_fitting(
  x_train_fun_array[train_train_ind,,,],
  y_train[train_train_ind,],
  list(x_train_fun_array[train_val_ind,,,], 
       y_train[train_val_ind,])
)

prediction_imagenet <- fit_imagenet %>% predict(x_test_fun_array)

saveRDS(prediction_imagenet, file="results/prediction_imagenet.RDS")

rm(fit_imagenet, prediction_imagenet); gc()

################ Convolutional Neural Network ##################

fit_cnn <- conv_fitting(
  x_train_fun_array[train_train_ind,,,],
  y_train[train_train_ind,],
  list(x_train_fun_array[train_val_ind,,,], 
       y_train[train_val_ind,])
)

prediction_cnn <- fit_cnn %>% predict(x_test_fun_array)

saveRDS(prediction_cnn, file="results/prediction_cnn.RDS")

rm(fit_cnn, prediction_cnn); gc()

################ GPFDA ####################
# vignette("gpfr", package = "GPFDA")
# package does not work with the additional array class for matrices
xtf <- lapply(x_train_fun, function(x){ class(x) <- "matrix"; return(x)})

fit_gp <- gpfr(response = y_train, 
               time = 1:101, 
               uReg = x_train_sca[,22:25],
               fxReg = x_train_fun,
               gpReg = matrix(rep(1:101, nrow(y_train)), ncol=101),
               fyList = list(nbasis = 15, lambda = 0.0001),
               fxList = list(list(nbasis = 15, lambda = 0.0001))[rep(1,5)],
               uCoefList = list(list(nbasi = 6, lambda = 0.0001)),
               Cov = 'pow.ex', gamma = 1, fitting = T)

prediction_gp <- predict(object=fit_gp, 
                         fxReg = x_test_fun, 
                         testInputGP = matrix(rep(1:101, nrow(y_test)), ncol=101),
                         testTime = 1:101,
                         uReg = x_test_sca[,22:25], 
                         gpReg = NULL
)


saveRDS(prediction_gp, file="results/prediction_gp.RDS")

rm(fit_gp, prediction_gp); gc()

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

train <- as.data.frame(lapply(train, function(x) if(!is.null(dim(x))) return(I(x)) else x))
train$id <- as.factor(train$id)
train$cond <- as.factor(train$cond)
train$sex <- as.factor(train$sex)
test <- as.data.frame(lapply(test, function(x) if(!is.null(dim(x))) return(I(x)) else x))
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

saveRDS(prediction_fdboost, file="results/prediction_fdboost.RDS")

rm(m, prediction_fdboost); gc()

################ FunFor ####################

datatrain <- cbind(do.call("cbind", x_train_fun), x_train_sca)
colnames(datatrain) <- paste0("X", 1:ncol(datatrain))
datatest <- cbind(do.call("cbind", x_test_fun), x_test_sca)
colnames(datatest) <- paste0("X", 1:ncol(datatest))

# apply minimal jitter because FunFor will otherwise fail
datatrain <- datatrain + matrix(rnorm(prod(dim(datatrain)), mean=0, sd=1e-8), ncol=ncol(datatrain))
datatest <- datatest + matrix(rnorm(prod(dim(datatest)), mean=0, sd=1e-8), ncol=ncol(datatest))

colnames(y_train) <- paste0("Y", 1:101)

datatrain = cbind(as.data.frame(y_train), as.data.frame(datatrain))

form = paste("datatrain[, 1:101]", "~", paste(colnames(datatrain[,-1*(1:101)]), collapse = " + "))

# o_split = FunFor:::optimal_size(form, datatrain) # too costly
funfor_fit = FunFor:::FunFor(form, datatrain, mtry = 40, ntree = 10, npc = 5, m_split = 8)

prediction_funfor <- FunFor:::predict.mvRF(funfor_fit, as.data.frame(datatest))

saveRDS(prediction_funfor, file="results/prediction_funfor.RDS")

rm(funfor_fit, prediction_funfor); gc()

################ Functional Intercept ####################

response <- response_vars[1]

# Create formula
form <- paste(response, " ~ 1")

# initialize the model
mint <- pffr(as.formula(form),
          yind = 1:101,
          algorithm = "bam",
          data = train)

prediction_intercept <- mint %>% predict(test)

saveRDS(prediction_intercept, file="results/prediction_intercept.RDS")

rm(mint, prediction_intercept); gc()

################ deepregression ####################

################ Implicit approach (XGBoost) ####################

# xgb_obj <- "reg:squarederror"
# 
# xgb_data <- xgb.DMatrix(data = as.matrix(data[rep(1:nrow(data[train_train_ind,]), each=101),]), 
#                         label = c(y_train[train_train_ind,]))
# 
# xgb_val <- xgb.DMatrix(data = as.matrix(data[rep(1:nrow(data[train_val_ind,]), each=101),]), 
#                         label = c(y_train[train_val_ind,]))
# 
# xgb_mod <- xgb.train(data = xgb_data, 
#                      early_stopping_rounds = 10, 
#                      watchlist = list(val=xgb_val),
#                      nrounds = 5000, 
#                      params = list(objective = xgb_obj,
#                                    max_depth = 30),
#                      verbose = FALSE)
# 
# prediction_xgb <- predict(xgb_mod, 
#                           newdata = xgb.DMatrix(data = as.matrix(data[rep(1:nrow(datatest), 
#                                                                           each=101),])))
# prediction_xgb <- matrix(prediction_xgb, ncol=101, byrow = F)
# 
# saveRDS(prediction_xgb, file="results/prediction_xgb.RDS")
# 
# rm(xgb_mod, prediction_xgb); gc()

###############################################################
######################### Comparison ##########################
###############################################################

results <- c(list(y_test), lapply(list.files("results", full.names = T), function(x) as.matrix(readRDS(x))))
nams <- c("truth", gsub("prediction_(.*)\\.RDS", "\\1", list.files("results")))
resultsDF <- do.call("rbind", lapply(1:length(results), function(i) 
  data.frame(value = c(results[[i]]), time = rep(1:101, each = nrow(y_test)),
             obs = rep(1:nrow(y_test), 101), what = nams[i])))

random_g <- mean(sapply(1:100, function(i) 
  mean(sqrt(apply((y_train[sample(1:nrow(y_train), nrow(y_test)),]-results[[1]])^2, 1, sum)))))

rmseDF <- data.frame(rmse = paste0("RMSE: ", round(
  sapply(2:length(results), function(i) mean(sqrt(apply((results[[i]]-results[[1]])^2, 1, sum)))), 4)),
                     what = nams[-1], obs = 1)

meth_lab <- c("CNN", "FDboost", "FuncNN", "FunFor", "ImageNet", "Intercept", "pffr")

rmseDF$what <- factor(rmseDF$what, levels = unique(rmseDF$what),
                      labels = meth_lab)
resultsDF$what <- factor(resultsDF$what, levels = unique(resultsDF$what),
                         labels = c("Truth", meth_lab))

ggplot(resultsDF %>% filter(!what %in% c("Intercept", "FunFor")), 
       aes(x = time, colour = what, y = value, group=obs)) + 
  geom_path(alpha=0.25) + theme_bw() + facet_wrap(~what) + ylim(-1.3,1.3) + 
  guides(colour = "none") + geom_text(
    data    = rmseDF %>% filter(!what %in% c("Intercept", "FunFor")),
    colour = "black",
    mapping = aes(x = -Inf, y = -Inf, label = rmse),
    hjust   = -0.1,
    vjust   = -1
  )

# ggsave(filename = "ml_comparison.pdf", width=7, height=5)

# Filter the results
filtered_results <- subset(resultsDF, !what %in% c("Intercept", "FunFor"))
filtered_rmse <- subset(rmseDF, !what %in% c("Intercept", "FunFor"))

# Plotting each method in a separate plot
palette <- RColorBrewer::brewer.pal(8, "Dark2")
palette <- col2rgb(palette)
palette <- rgb(palette[1,] / 255, palette[2,] / 255, palette[3,] / 255, alpha = 0.25)

pdf("ml_comparison.pdf", width = 10, height = 6)
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
dev.off()
