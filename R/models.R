# mistyR model training functions
# Copyleft (ɔ) 2020 Jovan Tanevski <jovan.tanevski@uni-heidelberg.de>


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
build_model <- function(views, target, model.function, model.name,
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
  model.views <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map(function(view) {
      
      model.view.cache.file <-
        paste0(
          cache.location, .Platform$file.sep,
          "model.", model.name, ".", view[["abbrev"]], ".", target,
          ".par.", cv.folds, ".", 
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
        
        # Build model
        model.view <- model.function(view_data = transformed.view.data,
                                     target = target, seed = seed, ...
                                     )
        
        if (cached) {
          readr::write_rds(model.view, model.view.cache.file)
        }
      }
      
      return(model.view)
    })

  jit <- withr::with_seed(
    seed,
    stats::rnorm(length(target.vector), 0, .Machine$double.eps)
  )

  # make oob predictions
  oob.predictions <- model.views %>%
    purrr::map(~ .x$unbiased.predictions$prediction) %>%
    rlist::list.cbind() %>%
    tibble::as_tibble(.name_repair = make.names) %>%
    dplyr::mutate(
      # add jitter term not only if sd(.x) == 0, but also if there are fewer 
      # unique values than cv.folds (more strict)
      # added this in case sd != 0, but there mostly 0, and say ten 1s, 
      # then chances are hight than one cv folds below will only contain 0s
      dplyr::across(where(~ length(unique(.x)) < cv.folds), ~ .x + jit),
      !!target := target.vector
    )
  
  # train lm on above, if bypass.intra set intercept to 0
  formula <- stats::as.formula(
    ifelse(bypass.intra, paste0(target, " ~ 0 + ."), paste0(target, " ~ ."))
  )
  
  # added 2. line to prevent "Error in svd(X) : infinite or missing values in 'x'"
  # occuring in ridge::linearRidge if one column has only 1 unique value
  #if (ncol(oob.predictions) <= 2 | 
  #    sum(map_lgl(oob.predictions, ~ length(unique(.x)) < 2)) > 0) {
  if (ncol(oob.predictions) <= 2) {
    combined.views <- stats::lm(
      formula,
      oob.predictions
    )
  } else {
    combined.views <- ridge::linearRidge(formula,
       oob.predictions,
       lambda = "automatic",
       scaling = "corrForm"
    )
  }

  # cv performance estimate
  test.folds <- withr::with_seed(
    seed,
    caret::createFolds(target.vector, k = cv.folds)
  )
  
  intra.view.only <-
    model.views[["intraview"]]$unbiased.predictions$prediction %>%
    tibble::enframe(name = NULL, value = "intraview") %>%
    dplyr::mutate(!!target := target.vector)
  
  performance.estimate <- test.folds %>% purrr::map_dfr(function(test.fold) {
    meta.intra <- stats::lm(
      formula,
      intra.view.only %>% dplyr::slice(-test.fold),
    )
    
    if (identical(oob.predictions, intra.view.only)) {
      meta.multi <- meta.intra
    } else {
      
      meta.multi <- ridge::linearRidge(
        formula,
        oob.predictions %>% dplyr::slice(-test.fold),
        lambda = "automatic",
        scaling = "corrForm"
      )
    }
    
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
      intra.RMSE = intra.RMSE, intra.R2 = 100 * intra.R2,
      multi.RMSE = multi.RMSE, multi.R2 = 100 * multi.R2
    )
  })
  
  final.model <- list(
    meta.model = combined.views,
    model.importances = purrr::map(model.views, ~ .x$importances),
    performance.estimate = performance.estimate
  )
  
  return(final.model)
}
