# mistyR model training functions
# Copyleft (ɔ) 2020 Jovan Tanevski <jovan.tanevski@uni-heidelberg.de>

#' Generate indices for bootstrap samples (bags)
#'
#' Helper function
#'
#' @param input expression matrix, rows correspond to observations (spots, cells),
#' columns correspond to features (genes).
#' @param n.bags number of bootstrap samples to generate.
#' @param seed seed
#'
#' @return list of length \item{\var{n.bags}} with indices for in bag
#' and out of bag (OOB) samples
#'
#' @noRd
generate_bags <- function(input, n.bags, seed) {
  seq.input <- seq_len(nrow(input))
  raw.bags <- withr::with_seed(
    seed,
    caret::createResample(seq.input, times = n.bags)
  )
  
  purrr::map(raw.bags, function(bag) {
    in.bag <- bag
    out.bag <- dplyr::setdiff(seq.input, in.bag)
    list(x.in = in.bag, x.out = out.bag)
  })
}

#' Generate indices for k-fold cross validation (CV)
#'
#' Helper function
#'
#' @param input expression matrix, rows correspond to observations (spots, cells),
#' columns correspond to features (genes).
#' @param n.folds k
#' @param seed seed
#'
#' @return list of length \item{\var{n.bags}} with indices for within fold
#' and holdout observations.
#'
#' @noRd
generate_folds <- function(input, n.folds, seed) {
  
  assertthat::assert_that(nrow(input) >= n.folds,
    msg = "The data has less rows than the requested number of cv folds."
  )
  
  seq.input <- seq_len(nrow(input))
  raw.folds <- withr::with_seed(
    seed,
    caret::createFolds(seq.input, k = n.folds)
  )
  
  purrr::map(raw.folds, function(fold) {
    hold.out <- fold
    in.fold <- dplyr::setdiff(seq.input, hold.out)
    list(x.in = in.fold, x.out = hold.out)
  })
}

#' Training a single model
#'
#' Helper function to train one (base) model, except for when random forest 
#' (ranger) is the learner, then ranger is directly called (generation of
#' bags and training of decision trees) is handled internally
#'
#' @param input expression matrix, rows correspond to observations (spots, cells),
#' columns correspond to features (genes).
#' @param target target, which must correspond to a column of the input.
#' @param learner which ML model to use.
#' @param n.bags how many bags to generate (and thus number of base learners).
#' Only on this function to be able to call ranger with the right number
#' of trees.
#' @param subs indices of one bag or CV folds.
#' @param n.vars how many variable to use for prediction (feature subspace
#' selection). Only in this function to be able to call ranger with the 
#' right number for mtry.
#' @param vars which variable to use for prediction.
#' @param seed seed.
#' @param ... feed additional parameters for the chosen ML model, e.g. for 
#' ranger splitrule = "extratrees".
#'
#' @return list containing the ML model (type of the object depending 
#' on which model was used) and unbiased predictions (for bags, the out-of-bag
#' predictions, and for CV, the prediction for the holdout set).
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
           assertthat::assert_that(requireNamespace("kernlab", quietly = TRUE),
            msg = "The package kernlab is required to use linear SVM"
           )
           
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
           assertthat::assert_that(requireNamespace("earth", quietly = TRUE),
            msg = "The package earth is required to use MARS."
           )
           
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

#' Training a bagged model
#'
#' Helper function to train a bagged model
#'
#' @param input expression matrix, rows correspond to observations (spots, cells),
#' columns correspond to features (genes).
#' @param target target, which must correspond to a column of the input.
#' @param learner which ML model to use.
#' @param n.bags how many bags to generate (and thus number of base learners).
#' Only on this function to be able to call ranger with the right number
#' of trees.
#' @param n.vars how many variable to use for prediction (feature subspace
#' selection). Only in this function to be able to call ranger with the 
#' right number for mtry.
#' @param seed seed.
#' @param ... feed additional parameters for the chosen ML model, e.g. for 
#' ranger splitrule = "extratrees".
#'
#' @return Depending on the learner list contraining: 
#' 1a) a single model (if \item{\var{learner}} == "ranger") or 
#' 1b) a list of length \item{\var{n.bags}} of models, and
#' 2) a tibble containing the unbiased predictions.
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

#' Training a model based on cross validation (CV)
#'
#' Helper function to train a model based on k-fold cross validation
#'
#' @param input expression matrix, rows correspond to observations (spots, cells),
#' columns correspond to features (genes).
#' @param target target, which must correspond to a column of the input.
#' @param learner which ML model to use.
#' @param cv.folds k 
#' @param seed seed.
#' @param ... feed additional parameters for the chosen ML model, e.g. for 
#' ranger splitrule = "extratrees".
#'
#' @return list containing the ML model (type of the object depending 
#' on which model was used) and unbiased predictions (for bags, the out-of-bag
#' predictions, and for CV, the prediction for the holdout set).
#'
#' @noRd
build_cv_model <- function(input, target, learner, cv.folds, seed, ...) {
  
  set.seed(seed)
  
  cvs <- generate_folds(input = input, n.folds = cv.folds, seed = seed)
  vars <- colnames(input)[colnames(input) != target]
  
  models <- purrr::map(cvs, function(cv) {
    
    model_wrapper(input = input, target = target, learner = learner, 
                  n.bags = n.bags, subs = cv, vars = vars, n.vars = NULL, 
                  seed = seed, ... = ...)
  })
  
  predictions <- purrr::map2_dfr(cvs, models, function(fold, model) {
    tibble::tibble(index = fold$x.out, prediction = model$prediction)
  }) %>% dplyr::arrange(index)
  
  w.model <- model_wrapper(input, target, learner = learner, n.bags = n.bags, 
                           n.vars = n.vars, seed = seed, 
                           subs = list(x.in = seq_len(nrow(input)), 
                                       x.out = seq_len(nrow(input))),
                           vars = vars, ... = ...)$model
  
  list(unbiased.predictions = predictions, 
       models = list(model = w.model))
}

#' Merging two named lists
#' 
#' Helper function so we do not have to rely on rlist::list.merge which
#' removes all list entries whose value is NULL
#' 
#' @param l1 list 1, whose values are potentially overwritten by \item{\var{l2}}
#' @param l2 list 2
#' 
#' @noRd
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

#' Helper function to extract importances from different models
#'
#' Since different ML algorithms can be used to model the different views, 
#' this function is needed to extract the importances from the model used.
#' 
#' @param models list containing the model (for ranger and CV based models) 
#' or models (for bagged models)
#' @param learner string corresponding to the ML models which ML model was used
#'
#' @return named vector containing the unnormalized importances for each
#' predictor
#'
#' @noRd
imp_model <- function(models, learner) {
  
  switch(learner,
         "ranger" = {
           models[[1]]$variable.importance
         },
         "lm" = {
           if (length(models) == 1) {
             coefs <- models[[1]]$coefficients
             coefs[names(coefs) != "(Intercept)"]
           } else {
             coefs <- purrr::map_dfr(models, function(model) {
               model$coefficients }) %>%
               colMeans(na.rm = TRUE)
             coefs[names(coefs) != "(Intercept)"]
           }
         },
         "svmLinear" = {
           if (length(models) == 1) {
             (t(models[[1]]@coef) %*% models[[1]]@xmatrix)[1,]
           } else {
             purrr::map_dfr(models, function(model) {
               # scaling or no scaling?
               # t(m@coef) %*% as.matrix(expr[bag$in.bag, -1][m@SVindex,])
               coefs <- t(model@coef) %*% model@xmatrix
               names(coefs) <- colnames(coefs)
               coefs
             } ) %>% colMeans(na.rm = TRUE)
           }
         },
         "earth" = {
           if (length(models) == 1) {
             coefs <- earth::evimp(models[[1]], trim = FALSE, sqrt. = TRUE)[, 6]
             names(coefs) <- stringr::str_remove(names(coefs), "-unused")
             coefs
           } else {
             purrr::map_dfr(models, function(model) {
               coefs <- earth::evimp(model, trim = FALSE, sqrt. = TRUE)[, 6] 
               names(coefs) <- stringr::str_remove(names(coefs), "-unused")
               coefs
             }) %>% colMeans(na.rm = TRUE)
           }
         }
  ) 
}

#' Train a multi-view model for a single target
#'
#' Trains individual models for each views. The models can either be 
#' ensemble models based on bagging (\item{\var{method}} = "bag") with the following base learners: 
#' decision trees ("ranger"), linear model ("lm"), support vector
#' machine with linear kernel ("linearSVM"), or multivariate adaptive
#' regression splines ("earth"). Or the models can be based on k-fold 
#' cross validation (\item{\var{method}} = "cv"), at the moment the same
#' models as for bagging are implemented.
#' 
#' The unbiased predictions from the view-specfic models, meaning the out-of-bag
#' or holdout predictions are then used to train linear meta model whose
#' overall performance is estimated by cross-validation.
#'
#' Default model is \code{\link[ranger]{ranger}()} for training the
#' view-specific models with the following parameters: \code{num.trees = 100}, 
#' \code{importance = "impurity"}, \code{num.threads = 1}, \code{seed = seed}.
#'
#' @inheritParams run_misty
#'
#' @param target name of the target marker to train models for.
#'
#' @return A list containing the trained meta-model, the view-specific 
#' importances and performance estimates.
#'
#' @noRd
build_model <- function(views, target, method, learner, n.vars, n.learners, 
                        cv.folds, bypass.intra, seed, cached, ...) {
  
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
  # so we are actually returning a model for each view by mapping over the
  # views
  model.views <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map(function(view) {
      model.view.cache.file <-
        paste0(
          cache.location, .Platform$file.sep,
          "model.", method, ".", learner, ".", view[["abbrev"]], ".", target,
          ".par", n.learners, ".", cv.folds, ".", n.vars, ".", 
          ellipsis.args.text, ".rds"
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
        
          model.view <- build_bagged_model(
            input = transformed.view.data, target = target, learner = learner, 
            n.bags = n.learners, n.vars = n.vars, seed = seed, ... = ...
            )
          
        } else if (method == "cv") {
          
          model.view <- build_cv_model(
            input = transformed.view.data, target = target, learner = learner,
            cv.folds = cv.folds, seed = seed, ... = ...
            )
        }
        
        if (cached) {
          readr::write_rds(model.view, model.view.cache.file)
        }
      }
      
      return(model.view)
    })

  # get oob predictions
  oob.predictions <- model.views %>%
    purrr::map_dfc(~ .x$unbiased.predictions$prediction) %>%
    dplyr::mutate(!!target := target.vector)

  # train lm on above, if bypass.intra set intercept to 0
  formula <- stats::as.formula(
    ifelse(bypass.intra, paste0(target, " ~ 0 + ."), paste0(target, " ~ ."))
  )
  combined.views <- stats::lm(
    formula,
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
      formula,
      intra.view.only %>% dplyr::slice(-test.fold)
    )
    meta.multi <- stats::lm(
      formula,
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
    # necessary? if target and prediction are all 0 we get NaN otherwise.
    intra.R2 <- ifelse(intra.RMSE == 0, 1, intra.R2)
    
    multi.RMSE <- caret::RMSE(multi.view.prediction, target.vector[test.fold])
    multi.R2 <- caret::R2(multi.view.prediction, target.vector[test.fold],
      formula = "traditional"
    )
    # necessary? if target and prediction are all 0 we get NaN otherwise.
    multi.R2 <- ifelse(intra.RMSE == 0, 1, multi.R2)
    
    tibble::tibble(
      intra.RMSE = intra.RMSE, intra.R2 = 100*intra.R2,
      multi.RMSE = multi.RMSE, multi.R2 = 100*multi.R2
    )
  })
  
  final.model <- list(
    meta.model = combined.views,
    model.importances = map(model.views, function(model.view) {
      imp_model(models = model.view$model, learner = learner) }),
    performance.estimate = performance.estimate
  )
  
  return(final.model)
}
