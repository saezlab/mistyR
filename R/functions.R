# MISTy view-specific functions
# Copyleft (ɔ) 2020 Jovan Tanevski [jovan.tanevski@uni-heidelberg.de]

#' function to merge named arguments from two lists
#' 
#' @export
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

#' RF implementation
#' 
#' @export
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

#' Bagged MARS
#' 
#' @export
bagged_earth_model <- function(view_data, target, seed, n.vars = NULL, 
                               n.learners = 100, ...) {
  
  ellipsis.args <- list(...)
  
  # get ellipsis arguments
  if ("n.vars" %in% ellipsis.args) n.vars <- ellipsis.args$n.vars
  
  if ("n.learners" %in% ellipsis.args) n.learners <- ellipsis.args$n.learners
  
  # how many predictors and how many variables to consider for each bag
  n.vars <- ifelse(is.null(n.vars), ncol(view_data)-1, n.vars)
  predictors <- colnames(view_data)[colnames(view_data) != target]
  
  # generate the bags
  bags <- withr::with_seed(
    seed,
    caret::createResample(1:nrow(view_data), times = n.learners)
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

#' Bagged Linear Model
#' 
#' @export
bagged_linear_model = function(view_data, target, seed, n.vars = NULL, 
                               n.learners = 100, ...) {
  
  ellipsis.args <- list(...)
  
  # get ellipsis arguments
  if ("n.vars" %in% ellipsis.args) n.vars <- ellipsis.args$n.vars
  
  if ("n.learners" %in% ellipsis.args) n.learners <- ellipsis.args$n.learners
  
  # how many predictors and how many variables to consider for each bag
  n.vars <- ifelse(is.null(n.vars), ncol(view_data)-1, n.vars)
  predictors <- colnames(view_data)[colnames(view_data) != target]
  
  # generate the bags
  bags <- withr::with_seed(
    seed,
    caret::createResample(1:nrow(view_data), times = n.learners)
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

#' Simple Linear Model
#' 
#' @export
linear_model = function(view_data, target, seed, k = 10, ...) {
  
  ellipsis.args <- list(...)
  
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


#' SVM Implementation
#' 
#' @export
svm_model = function(view_data, target, seed, approx = TRUE, 
                        approx.frac = 0.2, k = 10, ...) {
  
  ellipsis.args <- list(...)
  
  # get these parameters from the ellipsis
  if ("approx.frac" %in% ellipsis.args) approx.frac <- ellipsis.args$approx.frac
  if ("frac" %in% ellipsis.args) frac <- ellipsis.args$frac
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = k)
  )
  
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    # subsampling to reduce the computational cost
    if (approx) in.fold <- in.fold[sample(1:length(in.fold), 
                                          length(in.fold)*approx.frac)]
    
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


#' Boosted Trees Implementation
#' 
#' @export
boosted_trees_model = function(view_data, target, seed, k = 10, ...) {
  
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
      nrounds = 10,
      verbose = 0,
      booster = "gbtree", 
      eta = 0.3, 
      objective = "reg:squarederror",
      max_depth = 6, 
      nthread = 1
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_2(algo.arguments, ellipsis.args)
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
    nrounds = 10,
    verbose = 0, 
    booster = "gbtree", 
    eta = 0.3, 
    objective = "reg:squarederror",
    max_depth = 6, 
    nthread = 1
  )
  
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_2(algo.arguments.wm, ellipsis.args)
  }
  
  whole.model <- do.call(xgboost::xgboost, algo.arguments.wm)
  
  importance_matrix <- xgboost::xgb.importance(model = whole.model)
  importances <- unlist(importance_matrix[, "Gain"])
  names(importances) <- unlist(importance_matrix[, "Feature"])
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}


#' Neural Network Implementation
#' 
#' @export
nn_model = function(view_data, target, seed, k = 10, ...) {
  
  ellipsis.args <- list(...)
  
  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = k)
  )
  
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]
    
    train <- view_data[in.fold, ]
    test <- view_data[holdout, ]
    
    algo.arguments <- list(
      as.formula(paste0(target, "~ .")),
      data = train,
      hidden = c(10),
      linear.output = FALSE,
      lifesign = "none",
      rep = 1,
      stepmax = 1e4
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_2(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(neuralnet::neuralnet, algo.arguments)
    
    label.hat <- predict(model, test)
    
    tibble::tibble(index = holdout, prediction = label.hat)
  }) %>% dplyr::arrange(index)
  
  algo.arguments.wm <- list(
    as.formula(paste0(target, "~ .")),
    data = view_data,
    hidden = c(10),
    linear.output = FALSE,
    lifesign = "none",
    rep = 1,
    stepmax = 1e4
  )
  
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_2(algo.arguments.wm, ellipsis.args)
  }
  
  whole.model <- do.call(neuralnet::neuralnet, algo.arguments.wm)
  
  predictor <- Predictor$new(whole.model, 
                             data = view_data %>% dplyr::select(-tidyselect::all_of(target)) , 
                             y = view_data %>% dplyr::pull(tidyselect::all_of(target)))
  
  imp <- FeatureImp$new(predictor, loss = "mse")$results
  importances <- imp$importance
  names(importances) <- imp$feature
  
  list(unbiased.predictions = holdout.predictions, 
       importances = importances)
}




