# mistyR model training functions
# Copyright (c) 2020 Jovan Tanevski <jovan.tanevski@uni-heidelberg.de>

#' Helper function to generate bootstrap aggregation samples
#'
#' @return Row indices for the bags
#'
#' @noRd
generate_bags <- function(input, n.bags, seed) {
  set.seed(seed)
  seq_len(n.bags) %>%
    purrr::map(function(bag) {
      in.bag <- sample.int(nrow(input), nrow(input), replace=TRUE)
      out.bag <- dplyr::setdiff(seq_len(nrow(input)), in.bag)
      list(x.in = in.bag, x.out = out.bag)
    })
}

#' Helper function to generate folds for k-fold cross validation
#'
#' @return Row indices for the bags
#'
#' @noRd
generate_folds <- function(input, n.folds, seed) {
  
  n.rows <- nrow(input)
  
  assertthat::assert_that(n.rows >= n.folds,
    msg = "The data has less rows than the requested number of cv folds."
  )
  
  set.seed(seed)
  shuffled <- sample(seq_len(n.rows))
  steps <- c(ceiling(n.rows/n.folds), floor(n.rows/n.folds))
  corr <- n.folds * steps[1] - n.rows
  
  seq_len(n.folds) %>%
    purrr::map(function(fold) {
      
      if (fold <= corr) { index <- seq.int(
        (fold-1)*steps[1]+1, 
        fold*steps[1]
      ) 
      } else { index <- seq.int(
        corr*steps[1] + (fold-corr-1)*steps[2] + 1,
        corr*steps[1] + (fold-corr)*steps[2]
      ) }
      
      hold.out <- shuffled[index]
      in.fold <- shuffled[!(shuffled %in% hold.out)]
      
      list(x.in = in.fold, x.out = hold.out)
    })
}

#' Helper function to build the ML models
#'
#' @return model
#'
#' @noRd
model_wrapper <- function(input, target, learner, n.bags, subs,
                          n.vars, vars, seed, ...) {
  ellipsis.args <- list(...)
  switch(learner,
         "ranger" = {
           
           algo.arguments <- list(
             formula = stats::as.formula(paste0(target, " ~ .")),
             data = input,
             num.trees = n.bags,
             importance = "impurity",
             mtry = n.vars, 
             verbose = FALSE, 
             num.threads = 1,
             seed = seed)
           
           if (!(length(ellipsis.args) == 0)) {
             algo.arguments <- merge_2(algo.arguments, ellipsis.args)
           }
           
           model <- do.call(ranger::ranger, algo.arguments)
           
           predictions <- tibble::tibble(index = seq_len(nrow(input)), 
                                         prediction = model$predictions)
           
           list(unbiased.predictions = predictions, 
                models = list(model = model))
           
         }, "lm" = {
           
           algo.arguments <- list(
             formula = stats::as.formula(paste0(target, " ~ .")),
             data = input[subs$x.in, c(target, vars)]
           )
           
           if (!(length(ellipsis.args) == 0)) {
             algo.arguments <- merge_2(algo.arguments, ellipsis.args)
           }
           
           model <- do.call(stats::lm, algo.arguments)
           
           pred <- predict.lm(model, input[subs$x.out, vars])
           list(model = model, prediction = pred)
           
         }, "svmLinear" = {
           
           algo.arguments <- list(
             x = stats::as.formula(paste0(target, " ~ .")),
             data = input[subs$x.in, c(target, vars)],
             kernel = "vanilladot",
             C = 1,
             type = "eps-svr",
             kpar = list() # no hyperparameters for linear kernel
           )
           
           if (!(length(ellipsis.args) == 0)) {
             algo.arguments <- merge_2(algo.arguments, ellipsis.args)
           }
           
           model <- do.call(kernlab::ksvm, algo.arguments)
           
           pred <- kernlab::predict(model, input[subs$x.out, vars])
           list(model = model, prediction = pred)
           
         }, "earth" = {
           
           algo.arguments <- list(
             formula = stats::as.formula(paste0(target, " ~ .")),
             data = input[subs$x.in, c(target, vars)]
           )
           
           if (!(length(ellipsis.args) == 0)) {
             algo.arguments <- merge_2(algo.arguments, ellipsis.args)
           }
           
           model <- do.call(earth::earth, algo.arguments)
           
           pred <- predict(model, input[subs$x.out, vars])
           list(model = model, prediction = pred)
         }
  )
}

#' Building ML model based on Bagging (Ensemble Model)
#'
#' @return Relative Importances
#'
#' @noRd
build_bagged_model <- function(input, target, learner, n.bags, n.vars, 
                               seed, ...) {
  
  set.seed(seed)
  
  if (learner == "ranger") {
    
    return(
      model_wrapper(input = input, target = target, learner = learner, 
                    n.bags = n.bags, subs = NULL, n.vars = n.vars, 
                    vars = NULL, seed = seed, ...) 
    )
    
  } else {
    
    bags <- generate_bags(input = input, n.bags = n.bags, seed = seed)
    
    n.vars <- ifelse(is.null(n.vars), ncol(input)-1, n.vars)
    predictors <- colnames(input)[colnames(input) != target]
    
    models <- purrr::map(bags, function(bag) {
      vars <- sample(predictors, n.vars)
      
      model_wrapper(input = input, target = target, learner = learner, 
                    n.bags = n.bags, subs = bag, n.vars = n.vars, 
                    vars = vars, seed = seed, ...)
    })
    
    # Generate the OOB predictions, using the average for each sample
    predictions <- purrr::map2_dfr(bags, models, function(bag, model) {
      tibble::tibble(index = bag$x.out, prediction = model$prediction)
    }) %>% 
      dplyr::group_by(index) %>%
      dplyr::summarise(prediction = mean(prediction)) %>%
      dplyr::arrange(index)
    
    assertthat::assert_that(nrow(predictions) == nrow(input),
                            msg = "There are too few bags to get OOB predictions for all observations.
      Consider increasing the number of bags or using CV.")
    
    return(list(unbiased.predictions = predictions, 
                models = purrr::map(models, ~ .x$model)))
  }
}

#' Building a CV Model
#' 
#' TODO
#'
#' @return Relative Importances
#'
#' @noRd
build_cv_model <- function(input, target, learner, cv.folds, seed, ...) {
  
  set.seed(seed)
  
  cvs <- generate_folds(input, cv.folds, seed = seed)
  vars <- colnames(input)[colnames(input) != target]
  
  models <- purrr::map(cvs, function(cv) {
    
    model_wrapper(input = input, target = target, learner = learner, 
                  n.bags = n.bags, subs = cv, vars = vars, n.vars = NULL, 
                  seed = seed, ...)
  })
  
  predictions <- purrr::map2_dfr(cvs, models, function(fold, model) {
    tibble::tibble(index = fold$x.out, prediction = model$prediction)
  }) %>% dplyr::arrange(index)
  
  w.model <- model_wrapper(input, target, learner = learner, n.bags = n.bags, 
                           n.vars = n.vars, seed = seed, 
                           subs = list(x.in = seq_len(nrow(input)), 
                                       x.out = seq_len(nrow(input))),
                           vars = vars, ...)$model
  
  list(unbiased.predictions = predictions, 
       models = list(model = w.model))
}

#' Helper function to merge two names lists
#' 
#' Helper function so we do not have to rely on rlist::list.merge which
#' removes all list entries whose value is NULL
#'
#' @noRd
merge_2 <- function(l1, l2) {
  
  n1 <- names(l1)
  n2 <- names(l2)
  # only in n1
  diff <- n1[!(n1 %in% n2)]
  n1_list <- diff %>%
    purrr::set_names() %>%
    purrr::map(function(name) l1[[name]])
  
  # also in n2
  union <- n2[!(n2 %in% diff)]
  n2_list <- union %>%
    purrr::set_names() %>%
    purrr::map(function(name) l2[[name]])
  return(c(n1_list, n2_list))
}

#' Train a multi-view model for a single target
#'
#' Trains individual Random Forest models for each view, a linear meta-model
#' using the out-of-bag predictions of the view specific models and estimates
#' the overall performance by cross-validation.
#'
#' Default values passed to \code{\link[ranger]{ranger}()} for training the
#' view-specific models: \code{num.trees = 100}, \code{importance = "impurity"},
#' \code{num.threads = 1}, \code{seed = seed}.
#'
#' @inheritParams run_misty
#'
#' @param target name of the target marker to train models for.
#'
#' @return A list containing the trained meta-model, a list of trained
#' view-specific models and performance estimates.
#'
#' @noRd
build_model <- function(views, target, method = "bag", learner = "ranger", 
                        n.vars = NULL, n.learners = 100, cv.folds = 10,
                        bypass.intra = FALSE, seed = 42, cached = FALSE, ...) {
  
  cache.location <- R.utils::getAbsolutePath(paste0(
    ".misty.temp", .Platform$file.sep,
    views[["misty.uniqueid"]]
  ))
  
  if (cached && !dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  }
  
  expr <- views[["intraview"]][["data"]]
  
  target.vector <- expr %>% dplyr::pull(target)
  
  ellipsis.args <- list(...)
  ellipsis.args.text <- paste(names(ellipsis.args), ellipsis.args,
                              sep = ".", collapse = "."
  )
  
  # returns a list of models
  model.views <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map(function(view) {
      # Adjustment needed here (TODO: Caching)
      model.view.cache.file <-
        paste0(
          cache.location, .Platform$file.sep,
          "model.", view[["abbrev"]], ".", target,
          ".par", ellipsis.args.text, ".rds"
        )
      
      if (file.exists(model.view.cache.file) & cached) {
        model.view <- readr::read_rds(model.view.cache.file)
      } else {
        if ((view[["abbrev"]] == "intra") & bypass.intra) {
          transformed.view.data <-
            tibble::tibble(!!target := target.vector, ".novar" := 0)
        } else {
          transformed.view.data <- view[["data"]] %>%
            dplyr::mutate(!!target := target.vector)
        }
        
        if (method == "bag") {
          algo.arguments <- list(input = transformed.view.data, target = target, 
                                 learner = learner, n.bags = n.learners, 
                                 n.vars = n.vars, 
                                 seed = seed)
          
          if (!(length(ellipsis.args) == 0)) {
            algo.arguments <- merge_2(algo.arguments, ellipsis.args)
          }
          
          model.view <- do.call(build_bagged_model, algo.arguments)
          
        } else if (method == "cv") {
          algo.arguments <- list(input = transformed.view.data, target = target, 
                                 learner = learner, n.learners = n.learners,
                                 cv.folds = cv.folds, n.vars = n.vars, 
                                 seed = seed)
          if (!(length(ellipsis.args) == 0)) {
            algo.arguments <- merge_2(algo.arguments, ellipsis.args)
          }
          model.view <- do.call(build_cv_model, algo.arguments)
        }
        
        # TODO: Caching
        if (cached) {
          readr::write_rds(model.view, model.view.cache.file)
        }
      }
      
      return(model.view)
    })
  
  # make oob predictions
  oob.predictions <- model.views %>%
    purrr::map_dfc(~ .x$unbiased.predictions$prediction) %>%
    dplyr::mutate(!!target := target.vector)
  
  # train lm on above
  combined.views <- stats::lm(
    stats::as.formula(paste0(target, "~.")),
    oob.predictions
  )
  
  # cv performance estimate
  test.folds <- withr::with_seed(
    seed,
    caret::createFolds(target.vector, k = cv.folds)
  )
  
  intra.view.only <-
    model.views[["intraview"]]$unbiased.predictions$prediction %>%
    tibble::enframe(name = NULL) %>%
    dplyr::mutate(!!target := target.vector)
  
  performance.estimate <- test.folds %>% purrr::map_dfr(function(test.fold) {
    meta.intra <- stats::lm(
      stats::as.formula(paste0(target, "~.")),
      intra.view.only %>% dplyr::slice(-test.fold)
    )
    meta.multi <- stats::lm(
      stats::as.formula(paste0(target, "~.")),
      oob.predictions %>% dplyr::slice(-test.fold)
    )
    
    intra.prediction <- stats::predict(meta.intra, intra.view.only %>%
                                         dplyr::slice(test.fold))
    multi.view.prediction <- stats::predict(meta.multi, oob.predictions %>%
                                              dplyr::slice(test.fold))
    
    intra.RMSE <- caret::RMSE(intra.prediction, target.vector[test.fold])
    intra.R2 <- caret::R2(intra.prediction, target.vector[test.fold],
                          formula = "traditional"
    )
    
    multi.RMSE <- caret::RMSE(multi.view.prediction, target.vector[test.fold])
    multi.R2 <- caret::R2(multi.view.prediction, target.vector[test.fold],
                          formula = "traditional"
    )
    
    tibble::tibble(
      intra.RMSE = intra.RMSE, intra.R2 = intra.R2,
      multi.RMSE = multi.RMSE, multi.R2 = multi.R2
    )
  })
  
  final.model <- list(
    meta.model = combined.views,
    model.views = model.views,
    performance.estimate = performance.estimate
  )
  
  return(final.model)
}