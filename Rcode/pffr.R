
rm(list=ls())
gc()
# Load packages --------------------------------------------------------------
library(refund) # functional regression models for comparison
# data wrangling
library(tidyverse)

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

################ PFFR ####################

response <- response_vars[1]

# Create formula
form <- paste(response, " ~ 1 + ", paste(
  paste0("sff(", pred_vars,
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
