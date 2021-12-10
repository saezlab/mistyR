# MISTy view-specific functions
# Copyleft (ɔ) 2020 Jovan Tanevski [jovan.tanevski@uni-heidelberg.de]

#' function to merge named arguments from two lists
#' 
#' @export
merge_two <- function(l1, l2) {
  
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

#' RF implementation
#' 
#' @export
random_forest_model <- function(view_data, target, seed, ...) {
  
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
    algo.arguments <- merge_two(algo.arguments, ellipsis.args)
  }
  
  model <- do.call(ranger::ranger, algo.arguments)
  
  predictions <- tibble::tibble(index = seq_len(nrow(view_data)), 
                                prediction = model$predictions)
  
  list(unbiased.predictions = predictions, 
       importances = model$variable.importance)
}

#' Bagged MARS
#' 
#' @export
bagged_mars_model <- function(view_data, target, seed, 
                               n.bags = 50, ...) {
  
  assertthat::assert_that(requireNamespace("earth", quietly = TRUE),
    msg = "The package earth is required to use mars"
  )
  
  ellipsis.args <- list(...)

  # generate the bags
  bags <- withr::with_seed(
    seed,
    caret::createResample(1:nrow(view_data), times = n.bags)
  )
  
  # build one model for each bag, return oob predictions and importances
  models <- purrr::map(bags, function(bag) {
    
    algo.arguments <- list(
      formula = stats::as.formula(paste0(target, " ~ .")),
      data = view_data[bag, ],
      degree = 2
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_two(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(earth::earth, algo.arguments)
    
    oob <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% bag)]
    
    pred <- predict(model, view_data[oob, ])
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
    # fix for bypass intra
    names(coefs) <- ifelse(is.na(names(coefs)), "no.var", names(coefs))
    coefs
  }) %>% colMeans(na.rm = TRUE)
  
  
  
  list(unbiased.predictions = predictions, 
       importances = importances)
}

#' Bagged MARS
#' 
#' @export
mars_model = function(view_data, target, seed, approx = 1.0, k = 10, ...) {
  
  assertthat::assert_that(requireNamespace("earth", quietly = TRUE),
    msg = "The package earth is required to use mars"
  )
  
  ellipsis.args <- list(...)
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = k)
  )
  
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    # subsampling to reduce the computational cost
    if (approx != 1) in.fold <- in.fold[sample(1:length(in.fold), length(in.fold)*approx)]
    
    train <- view_data[in.fold, ]
    test <- view_data[holdout, ]
    
    algo.arguments <- list(
      formula = stats::as.formula(paste0(target, " ~ .")),
      data = train,
      degree = 2
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_two(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(earth::earth, algo.arguments)
    
    label.hat <- predict(model, test)
    
    tibble::tibble(index = holdout, prediction = label.hat)
  }) %>% dplyr::arrange(index)
  
  algo.arguments.wm <- list(
    formula = stats::as.formula(paste0(target, " ~ .")),
    data = view_data,
    degree = 2
  )
  
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_two(algo.arguments.wm, ellipsis.args)
  }
  
  whole.model <- do.call(earth::earth, algo.arguments.wm)
  
  importances <- earth::evimp(whole.model, trim = FALSE, sqrt. = TRUE)[, 6]
  names(importances) <- stringr::str_remove(names(importances), "-unused")
  # fix for bypass intra
  if (is.na(names(importances))) importances <- c("no.var" = 0)
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}

#' Simple Linear Model
#' 
#' @export
linear_model = function(view_data, target, seed, k = 10, ...) {
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = k)
  )
  
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    algo.arguments <- list(
      formula = stats::as.formula(paste0(target, " ~ .")),
      data = view_data[in.fold, ]
    )
    
    model <- do.call(stats::lm, algo.arguments)
    
    pred <- predict.lm(model, view_data[holdout, ])
    
    tibble::tibble(index = holdout, prediction = as.vector(pred))
  }) %>% dplyr::arrange(index)
  
  algo.arguments.wm <- list(
    formula = stats::as.formula(paste0(target, " ~ .")),
    data = view_data
  )
  
  whole.model <- do.call(stats::lm, algo.arguments.wm)
  
  importances <- whole.model$coefficients[names(whole.model$coefficients) != "(Intercept)"]
  # fix for bypass intra (replace NA with 0 for consistent behavior)
  importances <- ifelse(is.na(importances), 0, importances)
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}


#' SVM Implementation
#' 
#' @export
svm_model = function(view_data, target, seed, approx = 0.4, k = 10, ...) {
  
  assertthat::assert_that(requireNamespace("kernlab", quietly = TRUE),
    msg = "The package kernlab is required to use linear SVM"
  )
  
  ellipsis.args <- list(...)
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = k)
  )
  
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    # subsampling to reduce the computational cost
    if (approx != 1) in.fold <- in.fold[sample(1:length(in.fold), length(in.fold)*approx)]
    
    algo.arguments <- list(
      x = stats::as.formula(paste0(target, " ~ .")),
      data = view_data[in.fold, ],
      kernel = "vanilladot",
      C = 1,
      type = "eps-svr",
      kpar = list() # no hyperparameters for linear kernel
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_two(algo.arguments, ellipsis.args)
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
    algo.arguments.wm <- merge_two(algo.arguments.wm, ellipsis.args)
  }
  
  whole.model <- do.call(kernlab::ksvm, algo.arguments.wm)
  
  importances <- (t(whole.model@coef) %*% whole.model@xmatrix)[1,]
  # fix for bypass intra
  if (is.null(names(importances))) importances <- c("no.var" = 0)
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}


#' Linear Gradient Boosting Implementation
#' 
#' @export
gradient_boosting_model = function(view_data, target, seed, k = 10, ...) {
  
  assertthat::assert_that(requireNamespace("xgboost", quietly = TRUE),
    msg = "The package xgboost is required to use gradient boosting"
  )
  
  ellipsis.args <- list(...)
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = k)
  )
  
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    train <- view_data[in.fold, ]
    test <- view_data[holdout, ]
    
    pred.train <- train %>% dplyr::select(-tidyselect::all_of(target)) %>% as.matrix
    label.train <- train %>% dplyr::pull(tidyselect::all_of(target))
    
    algo.arguments <- list(
      data = pred.train,
      label = label.train,
      booster = "gbtree",
      nrounds = 10,
      verbose = 0,
      eta = 0.3, 
      objective = "reg:squarederror",
      max_depth = 6, 
      nthread = 1
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_two(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(xgboost::xgboost, algo.arguments)
    
    pred.test <- test %>% dplyr::select(-tidyselect::all_of(target)) %>% as.matrix
    
    label.hat <- predict(model, pred.test)
    
    tibble::tibble(index = holdout, prediction = label.hat)
  }) %>% dplyr::arrange(index)
  
  predictors <- view_data %>% dplyr::select(-tidyselect::all_of(target)) %>% as.matrix
  labels <- view_data %>% dplyr::pull(tidyselect::all_of(target))
  
  algo.arguments.wm <- list(
    data = predictors,
    label = labels,
    booster = "gbtree",
    nrounds = 10,
    verbose = 0, 
    eta = 0.3, 
    objective = "reg:squarederror",
    max_depth = 6, 
    nthread = 1
  )
  
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_two(algo.arguments.wm, ellipsis.args)
  }
  
  whole.model <- do.call(xgboost::xgboost, algo.arguments.wm)
  
  # if bypass intra is true, we need to catch the error
  importances <- tryCatch({
    importance_matrix <- xgboost::xgb.importance(model = whole.model)
    importances <- unlist(importance_matrix[, "Gain"])
    names(importances) <- unlist(importance_matrix[, "Feature"])
    importances
  }, error = function(cond){
    importances <- rep(0, ncol(predictors))
    names(importances) <- colnames(predictors)
    importances
  })
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}

#' @export
mlp_model = function(view_data, target, seed, approx = 0.6, k = 10, ...) {
  
  assertthat::assert_that(requireNamespace("RSNNS", quietly = TRUE),
    msg = "The package RSNNS is required to use mlp"
  )
  
  ellipsis.args <- list(...)
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = k)
  )
  
  predictors <- colnames(view_data)[colnames(view_data) != target]
  X <- view_data %>% dplyr::select(tidyselect::all_of(predictors)) %>% as.matrix()
  Y <- view_data %>% dplyr::pull(target)
  
  # made this an imap to track the folds!
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    # subsampling to reduce the computational cost
    if (approx != 1) in.fold <- in.fold[sample(1:length(in.fold), length(in.fold)*approx)]
    
    X_train <- X[in.fold, ] %>% as.matrix
    Y_train <- Y[in.fold]
    
    X_test <- X[holdout, ] %>% as.matrix
    Y_test <- Y[holdout]
    
    algo.arguments <- list(
      x = X_train,
      y = Y_train,
      size = c(10)
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_two(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(RSNNS::mlp, algo.arguments)
    
    label.hat <- predict(model, X_test)
    
    tibble::tibble(index = holdout, prediction = label.hat)
  }) %>% dplyr::arrange(index)
  
  algo.arguments.wm <- list(
    x = X,
    y = Y,
    size = c(10)
  )
  
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_two(algo.arguments.wm, ellipsis.args)
  }
  
  whole.model <- do.call(RSNNS::mlp, algo.arguments.wm)
  
  # fix for bypass intra (or in general if there is only a single predictor)
  if (ncol(X) == 1) {importances <- c("no.var" = 0)
  } else {
    predictor <- iml::Predictor$new(
      model = whole.model, 
      data = as.data.frame(X), 
      y = Y)
    
    imp <- iml::FeatureImp$new(predictor, loss = "mse")$results
    importances <- imp$importance
    names(importances) <- imp$feature
  }
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}

