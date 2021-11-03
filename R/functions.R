
library(mistyR)
data("synthetic")
library(tidyverse)

# function to merge arguments
merge_2 <- function(l1, l2) {
  
  n1 <- names(l1)
  n2 <- names(l2)
  diff <- n1[!(n1 %in% n2)]
  n1_list <- diff %>%
    purrr::set_names() %>%
    purrr::map(function(name) l1[[name]])
  
  union <- n2[!(n2 %in% diff)]
  n2_list <- union %>%
    purrr::set_names() %>%
    purrr::map(function(name) l2[[name]])
  return(c(n1_list, n2_list))
}

#####

# keep the ellipsis to pass some parameters to ranger
ranger_model <- function(view_data, target, seed, ...) {
  
  ellipsis.args <- list(...)
  
  # default ranger arguments
  algo.arguments <- list(
    formula = stats::as.formula(paste0(target, " ~ .")),
    data = view_data,
    num.trees = 100,
    importance = "impurity",
    mtry = NULL, 
    verbose = FALSE, 
    num.threads = 1,
    seed = seed)
  
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments <- merge_2(algo.arguments, ellipsis.args)
  }
  
  model <- do.call(ranger::ranger, algo.arguments)
  
  predictions <- tibble::tibble(index = seq_len(nrow(view_data)), 
                                prediction = model$predictions)
  
  list(unbiased.predictions = predictions, 
       importances = model$variable.importance)
}

#####

expr <- synthetic$synthetic1 %>% select(-c(row, col, type))
test <- ranger_model(expr, "ECM", splitrule = "extratrees", 42)
test$unbiased.predictions
test$importances

view_data <- expr

#####

bagged_earth_model <- function(view_data, target, seed, n.vars...) {
  
  ellipsis.args <- list(...)
  
  # how many predictors and how many variables to consider for each bag
  n.vars <- ifelse(is.null(n.vars), ncol(view_data)-1, n.vars)
  predictors <- colnames(view_data)[colnames(view_data) != target]
  
  # generate the bags
  bags <- withr::with_seed(
    seed,
    caret::createResample(1:nrow(view_data), times = 100)
  )
  
  # build one model for each bag, return oob predictions and importances
  models <- purrr::map(bags, function(bag) {
    
    vars <- sample(predictors, n.vars)
    
    algo.arguments <- list(
      formula = stats::as.formula(paste0(target, " ~ .")),
      data = view_data[bag, c(target, vars)]
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_2(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(earth::earth, algo.arguments)
    
    oob <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% bag)]
    
    pred <- predict(model, view_data[oob, vars])
    list(model = model, 
         prediction = tibble::tibble(index = oob, prediction = as.vector(pred)))
  })
  
  # Generate the OOB predictions, using the average for each sample
  predictions <- purrr::map_dfr(models, function(model) {
    tibble::tibble(model$prediction)
  }) %>% 
    dplyr::group_by(index) %>%
    dplyr::summarise(prediction = mean(prediction)) %>%
    dplyr::arrange(index)
  
  assertthat::assert_that(nrow(predictions) == nrow(view_data),
      msg = "There are too few bags to get OOB predictions for all observations.
      Consider increasing the number of bags or using CV.")
  
  importances <- purrr::map_dfr(models, function(model) {
    coefs <- earth::evimp(model$model, trim = FALSE, sqrt. = TRUE)[, 6] 
    names(coefs) <- stringr::str_remove(names(coefs), "-unused")
    coefs
  }) %>% colMeans(na.rm = TRUE)
  
  
  list(unbiased.predictions = predictions, 
       importances = importances)
}

#####

expr <- synthetic$synthetic1 %>% select(-c(row, col, type))
test <- bagged_earth_model(expr, "ECM", 42, NULL)
test$unbiased.predictions
test$importances

#####

bagged_linear_model = function(view_data, target, seed, n.vars, ...) {
  
  ellipsis.args <- list(...)
  
  # how many predictors and how many variables to consider for each bag
  n.vars <- ifelse(is.null(n.vars), ncol(view_data)-1, n.vars)
  predictors <- colnames(view_data)[colnames(view_data) != target]
  
  # generate the bags
  bags <- withr::with_seed(
    seed,
    caret::createResample(1:nrow(view_data), times = 100)
  )
  
  # build one model for each bag, return oob predictions and importances
  models <- purrr::map(bags, function(bag) {
    
    vars <- sample(predictors, n.vars)
    
    algo.arguments <- list(
      formula = stats::as.formula(paste0(target, " ~ .")),
      data = view_data[bag, c(target, vars)]
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_2(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(stats::lm, algo.arguments)
    
    oob <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% bag)]
    
    pred <- predict.lm(model, view_data[oob, vars])

    list(model = model, 
         prediction = tibble::tibble(index = oob, prediction = as.vector(pred)))
  })
  
  # Generate the OOB predictions, using the average for each sample
  predictions <- purrr::map_dfr(models, function(model) {
    tibble::tibble(model$prediction)
  }) %>% 
    dplyr::group_by(index) %>%
    dplyr::summarise(prediction = mean(prediction)) %>%
    dplyr::arrange(index)
  
  assertthat::assert_that(nrow(predictions) == nrow(view_data),
      msg = "There are too few bags to get OOB predictions for all observations.
      Consider increasing the number of bags or using CV.")
  
  importances <- purrr::map_dfr(models, function(model) {
    coefs <- purrr::map_dfr(models, function(model) {
      model$model$coefficients }) %>%
      colMeans(na.rm = TRUE)
    coefs[names(coefs) != "(Intercept)"]
  }) %>% colMeans(na.rm = TRUE)
  
  
  list(unbiased.predictions = predictions, 
       importances = importances)
}

#####

expr <- synthetic$synthetic1 %>% select(-c(row, col, type))
test <- bagged_linear_model(expr, "ECM", 42, NULL)
test$unbiased.predictions
test$importances

#####

cv_linear_model = function(view_data, target, seed, ...) {

  ellipsis.args <- list(...)
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = 10)
  )
  
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    algo.arguments <- list(
      formula = stats::as.formula(paste0(target, " ~ .")),
      data = view_data[in.fold, ]
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_2(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(stats::lm, algo.arguments)
    
    pred <- predict.lm(model, view_data[holdout, ])
    
    tibble::tibble(index = holdout, prediction = as.vector(pred))
  }) %>% dplyr::arrange(index)
  
  algo.arguments.wm <- list(
    formula = stats::as.formula(paste0(target, " ~ .")),
    data = view_data
  )
  
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_2(algo.arguments.wm, ellipsis.args)
  }
  
  whole.model <- do.call(stats::lm, algo.arguments.wm)
  
  importances <- whole.model$coefficients[names(whole.model$coefficients) != "(Intercept)"]
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}

#####

expr <- synthetic$synthetic1 %>% select(-c(row, col, type))
test <- cv_linear_model(expr, "ECM", 42)
test$unbiased.predictions
test$importances

#####

cv_svm_model = function(view_data, target, seed, ...) {
  
  ellipsis.args <- list(...)
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = 10)
  )
  
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    algo.arguments <- list(
      x = stats::as.formula(paste0(target, " ~ .")),
      data = view_data[in.fold, ],
      kernel = "vanilladot",
      C = 1,
      type = "eps-svr",
      kpar = list() # no hyperparameters for linear kernel
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_2(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(kernlab::ksvm, algo.arguments)
    
    pred <- kernlab::predict(model, view_data[holdout, ])
    
    tibble::tibble(index = holdout, prediction = as.vector(pred))
  }) %>% dplyr::arrange(index)
  
  algo.arguments.wm <- list(
    x = stats::as.formula(paste0(target, " ~ .")),
    data = view_data,
    kernel = "vanilladot",
    C = 1,
    type = "eps-svr",
    kpar = list() # no hyperparameters for linear kernel
  )
  
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_2(algo.arguments.wm, ellipsis.args)
  }
  
  whole.model <- do.call(kernlab::ksvm, algo.arguments.wm)
  
  importances <- (t(whole.model@coef) %*% whole.model@xmatrix)[1,]
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}

#####

expr <- synthetic$synthetic1 %>% select(-c(row, col, type))
test <- cv_svm_model(expr, "ECM", 42)
test$unbiased.predictions
test$importances

#####

cv_boosted_trees_model = function(view_data, target, seed, ...) {
  
  ellipsis.args <- list(...)
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = 10)
  )
  
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    train <- expr[in.fold, ]
    test <- expr[holdout, ]
    
    pred.train <- train %>% select(-all_of(target)) %>% as.matrix
    label.train <- train %>% pull(all_of(target))
    
    algo.arguments <- list(
      params = list(booster = "gbtree", eta = 0.3, objective = "reg:squarederror",
                    max.depth = 6, nthread = 1),
      data = pred.train,
      label = label.train,
      nrounds = 10,
      verbose = 0
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_2(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(xgboost::xgboost, algo.arguments)

    pred.test <- test %>% select(-all_of(target)) %>% as.matrix
    
    label.hat <- predict(model, pred.test)
    
    tibble::tibble(index = holdout, prediction = label.hat)
  }) %>% dplyr::arrange(index)
  

  
  whole.model <- do.call(kernlab::ksvm, algo.arguments.wm)
  
  importances <- (t(whole.model@coef) %*% whole.model@xmatrix)[1,]
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}

#####

# Learning about xgboost

predictors <- expr %>% select(-ECM) %>% as.matrix
labels <- expr %>% pull(ECM)

test <- xgboost::xgboost(
  params = list(booster = "gbtree", eta = 0.3, objective = "reg:squarederror",
                max.depth = 2, nthread = 1),
  data = predictors,
  label = labels,
  nrounds = 10
)

pred <- predict(test, predictors)

importance_matrix <- xgboost::xgb.importance(model = test)
importance_matrix

importance <- unlist(importance_matrix[, "Gain"])
names(importance) <- unlist(importance_matrix[, "Feature"])
importance

#####

# Hyperparameter Tuning

folds <- withr::with_seed(
  seed,
  caret::createFolds(seq.int(1, nrow(view_data)), k = 10)
)

hyper.test <- imap_dfr(folds, function(holdout, i) {
  
  in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
  
  train <- expr[in.fold, ]
  test <- expr[holdout, ]
  
  pred.train <- train %>% select(-ECM) %>% as.matrix
  label.train <- train %>% pull(ECM)
  
  model <- xgboost::xgboost(
    params = list(booster = "gbtree", eta = 0.3, objective = "reg:squarederror",
                  max.depth = 6, nthread = 1),
    data = pred.train,
    label = label.train,
    nrounds = 10,
    verbose = 0
  )
  
  pred.test <- test %>% select(-ECM) %>% as.matrix
  label.test <- test %>% pull(ECM)
  label.hat <- predict(model, pred.test)
  
  RMSE <- caret::RMSE(label.hat, label.test)
  R2 <- caret::R2(label.hat, label.test)
  
  tibble::tibble(cv = i, r2 = R2, rmse = RMSE)
  
}) %>%
  rbind(., tibble::tibble(cv = "Summary", r2 = mean(.$r2), rmse = mean(.$rmse))) 

hyper.test


