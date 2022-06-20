# MISTy view-specific modeling functions
# Copyleft (ɔ) 2021 Philipp Schäfer, Jovan Tanevski <jovan.tanevski@uni-heidelberg.de>

#' Random Forest Implementation
#' 
#' The default function to model each view based on random forest as implemented 
#' in \code{\link[ranger]{ranger}()}.
#' 
#' The following parameters are the default configuration: \code{num.trees = 100},
#'  \code{importance = "impurity"}, \code{num.threads = 1}, \code{seed = seed}.
#' 
#' Can be explicitly called via 
#' \code{run_misty(views, model.function = random_forest_model)}
#' 
#' @param view_data \code{tibble} containing the expression of the markers for
#' each cell (rows = cells, cols = markers)
#' @param target \code{string} indicating the target marker, passed by run_misty()
#' @param seed \code{integer} passed by run_misty()
#' @param ... all additional parameters are passed to the chosen ML model
#'
#' @noRd
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
    seed = seed
  )

  if (!(length(ellipsis.args) == 0)) {
    algo.arguments <- merge_two(algo.arguments, ellipsis.args)
  }

  model <- do.call(ranger::ranger, algo.arguments)

  predictions <- tibble::tibble(
    index = seq_len(nrow(view_data)),
    prediction = model$predictions
  )

  list(
    unbiased.predictions = predictions,
    importances = model$variable.importance
  )
}

#' Gradient Boosting Implementation
#'
#' Function that can be passed to \code{run_misty()} via \code{run_misty(views, 
#' model.function = gradient_boosting_model)} to model each view using gradient 
#' boosting. The algorithm is implemented in the \code{xgboost} R package.
#' 
#' The following parameters are the default configuration: \code{booster = "gbtree"},
#' \code{rounds = 10}, \code{objective = "reg:squarederror"}. Set  \code{booster}
#' to \code{"gblinear"} for linear boosting.
#'
#' @param view_data \code{tibble} containing the expression of the markers for
#' each cell (rows = cells, cols = markers)
#' @param target \code{string} indicating the target marker, passed by run_misty()
#' @param seed \code{integer} passed by run_misty()
#' @param k \code{integer} indicating the folds used in cross validation to 
#' compute the unbiased estimates
#' @param ... all additional parameters are passed to the chosen ML model
#'
#' @noRd
gradient_boosting_model <- function(view_data, target, seed, k = 10, ...) {
  
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

    pred.train <- train %>%
      dplyr::select(-tidyselect::all_of(target)) %>%
      as.matrix()
    label.train <- train %>% dplyr::pull(tidyselect::all_of(target))

    algo.arguments <- list(
      data = pred.train,
      label = label.train,
      booster = "gbtree",
      nrounds = 10,
      verbose = 0,
      objective = "reg:squarederror",
      nthread = 1
    )

    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_two(algo.arguments, ellipsis.args)
    }

    model <- do.call(xgboost::xgboost, algo.arguments)

    pred.test <- test %>%
      dplyr::select(-tidyselect::all_of(target)) %>%
      as.matrix()

    label.hat <- predict(model, pred.test)

    tibble::tibble(index = holdout, prediction = label.hat)
  }) %>% dplyr::arrange(index)

  predictors <- view_data %>%
    dplyr::select(-tidyselect::all_of(target)) %>%
    as.matrix()
  labels <- view_data %>% dplyr::pull(tidyselect::all_of(target))

  algo.arguments.wm <- list(
    data = predictors,
    label = labels,
    booster = "gbtree",
    nrounds = 10,
    verbose = 0,
    objective = "reg:squarederror",
    nthread = 1
  )

  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_two(algo.arguments.wm, ellipsis.args)
  }

  whole.model <- do.call(xgboost::xgboost, algo.arguments.wm)

  # if bypass intra is true, we need to catch the error due to "no.var" = 0
  importances <- tryCatch(
    {
      importance_matrix <- xgboost::xgb.importance(model = whole.model)
      importances <- unlist(importance_matrix[, 2])
      names(importances) <- unlist(importance_matrix[, 1])
      importances
    },
    error = function(cond) {
      importances <- rep(0, ncol(predictors))
      names(importances) <- colnames(predictors)
      importances
    }
  )

  list(
    unbiased.predictions = holdout.predictions,
    importances = importances
  )
}

#' Bagged MARS
#'
#' Function that can be passed to \code{run_misty()} via \code{run_misty(views, 
#' model.function = bagged_mars_model)} to model each view using bagged MARS,
#' meaning multivariate adaptive spline regression models trained with 
#' bootstrap aggregation samples. The algorithm is implemented in the \code{earth} 
#' R package.
#' 
#' The following parameters are the default configuration: \code{degree = 2}.
#' Furthermore 50 base learners are used by default (see \code{n.bags}).
#' 
#' @param view_data \code{tibble} containing the expression of the markers
#' for each cell (rows = cells, cols = markers)
#' @param target \code{string} indicating the target marker, passed by run_misty()
#' @param seed \code{integer} passed by run_misty()
#' @param n.bags \code{integer} indicating number of bags
#' (bootstrap aggregation samples) used for bagging
#' @param ... all additional parameters are passed to the chosen ML model
#'
#' @noRd
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
      degree = 2,
      fast.k = 5,
      thresh = 0.01
    )

    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_two(algo.arguments, ellipsis.args)
    }

    model <- do.call(earth::earth, algo.arguments)

    oob <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% bag)]

    pred <- predict(model, view_data[oob, ])
    list(
      model = model,
      prediction = tibble::tibble(index = oob, prediction = as.vector(pred))
    )
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
      Consider increasing the number of bags or using CV."
  )

  importances <- purrr::map_dfr(models, function(model) {
    coefs <- earth::evimp(model$model, trim = FALSE, sqrt. = TRUE)[, 6]
    names(coefs) <- stringr::str_remove(names(coefs), "-unused")
    # fix for bypass intra
    names(coefs) <- ifelse(is.na(names(coefs)), "no.var", names(coefs))
    coefs
  }) %>% colMeans(na.rm = TRUE)



  list(
    unbiased.predictions = predictions,
    importances = importances
  )
}

#' Bagged MARS
#'
#' Function that can be passed to \code{run_misty()} via \code{run_misty(views, 
#' model.function = mars_model)} to model each view using a multivariate adaptive
#' spline regression model. The algorithm is implemented in the \code{earth} R package.
#' 
#' The following parameters are the default configuration: \code{degree = 2}
#' 
#' @param view_data \code{tibble} containing the expression of the markers for
#' each cell (rows = cells, cols = markers)
#' @param target \code{string} indicating the target marker, passed by run_misty()
#' @param seed \code{integer} passed by run_misty()
#' @param k \code{integer} indicating the folds used in cross validation to
#' compute the unbiased estimates
#' @param approx \code{double} indicating the fraction of the training data used
#' for training in each fold
#' @param ... all additional parameters are passed to the chosen ML model
#'
#' @noRd
mars_model <- function(view_data, target, seed, approx = 1.0, k = 10, ...) {
  
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
    if (approx != 1) in.fold <- sample(in.fold, length(in.fold) * approx)

    train <- view_data[in.fold, ]
    test <- view_data[holdout, ]

    algo.arguments <- list(
      formula = stats::as.formula(paste0(target, " ~ .")),
      data = train,
      degree = 2,
      fast.k = 5,   # reducing fast.k and
      thresh = 0.01 # increasing thresh to increase the performance
    )

    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_two(algo.arguments, ellipsis.args)
    }

    model <- do.call(earth::earth, algo.arguments)

    label.hat <- predict(model, test)

    tibble::tibble(index = holdout, prediction = label.hat[,1])
  }) %>% dplyr::arrange(index)

  algo.arguments.wm <- list(
    formula = stats::as.formula(paste0(target, " ~ .")),
    data = view_data,
    degree = 2,
    fast.k = 5,
    thresh = 0.01
  )

  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_two(algo.arguments.wm, ellipsis.args)
  }

  whole.model <- do.call(earth::earth, algo.arguments.wm)

  importances <- earth::evimp(whole.model, trim = FALSE, sqrt. = TRUE)[, 6]
  names(importances) <- stringr::str_remove(names(importances), "-unused")
  # fix for bypass intra
  if (all(is.na(names(importances)))) importances <- c("no.var" = 0)

  list(
    unbiased.predictions = holdout.predictions,
    importances = importances
  )
}

#' Simple Linear Model
#'
#' Function that can be passed to \code{run_misty()} via \code{run_misty(views, 
#' model.function = linear_model)} to model each view using a simple linear model
#' 
#' @param view_data \code{tibble} containing the expression of the markers for
#' each cell (rows = cells, cols = markers)
#' @param target \code{string} indicating the target marker, passed by run_misty()
#' @param seed \code{integer} passed by run_misty()
#' @param k \code{integer} indicating the folds used in cross validation to compute the
#' unbiased estimates
#' @param ... all additional parameters are passed to the chosen ML model
#'
#' @noRd
linear_model <- function(view_data, target, seed, k = 10, ...) {
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

  list(
    unbiased.predictions = holdout.predictions,
    importances = importances
  )
}


#' SVM Implementation
#'
#' Function that can be passed to \code{run_misty()} via \code{run_misty(views, 
#' model.function = gradient_boosting_model)} to model each view using a 
#' support vector machine. The algorithm is implemented in the \code{kernlab} R package.
#' 
#' The following parameters are the default configuration: \code{kernel = "vanilladot"} 
#' (linear kernel), \code{C = 1}, \code{type = "eps-svr"}.
#' 
#' @param view_data \code{tibble} containing the expression of the markers for each cell
#' (rows = cells, cols = markers)
#' @param target \code{string} indicating the target marker, passed by run_misty()
#' @param seed \code{integer} passed by run_misty()
#' @param k \code{integer} indicating the folds used in cross validation to compute the
#' unbiased estimates
#' @param ... all additional parameters are passed to the chosen ML model
#'
#' @noRd
svm_model <- function(view_data, target, seed, approx = 0.4, k = 10, ...) {
  
  assertthat::assert_that(requireNamespace("kernlab", quietly = TRUE),
    msg = "The package kernlab is required to use SVM"
  )

  ellipsis.args <- list(...)

  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = k)
  )

  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]

    # subsampling to reduce the computational cost
    if (approx != 1) in.fold <- sample(in.fold, length(in.fold) * approx)

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

  importances <- (t(whole.model@coef) %*% whole.model@xmatrix)[1, ]
  # fix for bypass intra
  if (is.null(names(importances))) importances <- c("no.var" = 0)

  list(
    unbiased.predictions = holdout.predictions,
    importances = importances
  )
}

#' Multi-layer Perceptron Implementation
#'
#' Function that can be passed to \code{run_misty()} via \code{run_misty(views, 
#' model.function = mlp_model)} to model each view using a multi-layer perceptron.
#'  The algorithm is implemented in the \code{RSNNS} R package.
#' 
#' The following parameters are the default configuration: \code{size = c(10)} 
#' (meaning we have 1 hidden layer with 10 units).
#' 
#' @param view_data \code{tibble} containing the expression of the markers for each cell
#' (rows = cells, cols = markers)
#' @param target \code{string} indicating the target marker, passed by run_misty()
#' @param seed \code{integer} passed by run_misty()
#' @param k \code{integer} indicating the folds used in cross validation to compute the
#' unbiased estimates
#' @param ... all additional parameters are passed to the chosen ML model
#'
#' @noRd
mlp_model <- function(view_data, target, seed, approx = 0.6, k = 10, ...) {
  
  assertthat::assert_that(requireNamespace("RSNNS", quietly = TRUE),
    msg = "The package RSNNS is required to use mlp"
  )

  ellipsis.args <- list(...)

  folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.int(1, nrow(view_data)), k = k)
  )

  predictors <- colnames(view_data)[colnames(view_data) != target]
  X <- view_data %>%
    dplyr::select(tidyselect::all_of(predictors)) %>%
    as.matrix()
  Y <- view_data %>% dplyr::pull(target)

  # made this an imap to track the folds!
  holdout.predictions <- purrr::map_dfr(folds, function(holdout) {
    in.fold <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% holdout)]

    # subsampling to reduce the computational cost
    if (approx != 1) in.fold <- sample(in.fold, length(in.fold) * approx)

    X_train <- X[in.fold, ] %>% as.matrix()
    Y_train <- Y[in.fold]

    X_test <- X[holdout, ] %>% as.matrix()
    Y_test <- Y[holdout]

    algo.arguments <- list(
      x = X_train,
      y = Y_train,
      size = c(10),
      learnFunc = "BackpropBatch"
    )

    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_two(algo.arguments, ellipsis.args)
    }

    model <- do.call(RSNNS::mlp, algo.arguments)

    label.hat <- predict(model, X_test)

    tibble::tibble(index = holdout, prediction = label.hat[,1])
  }) %>% dplyr::arrange(index)

  algo.arguments.wm <- list(
    x = X,
    y = Y,
    size = c(10),
    learnFunc = "BackpropBatch"
  )

  if (!(length(ellipsis.args) == 0)) {
    algo.arguments.wm <- merge_two(algo.arguments.wm, ellipsis.args)
  }

  whole.model <- do.call(RSNNS::mlp, algo.arguments.wm)

  # fix for bypass intra (or in general if there is only a single predictor)
  if (ncol(X) == 1) {
    importances <- c("no.var" = 0)
  } else {
    predictor <- iml::Predictor$new(
      model = whole.model,
      data = as.data.frame(X),
      y = Y
    )

    imp <- iml::FeatureImp$new(predictor, loss = "mse")$results
    importances <- imp$importance
    names(importances) <- imp$feature
  }

  list(
    unbiased.predictions = holdout.predictions,
    importances = importances
  )
}
