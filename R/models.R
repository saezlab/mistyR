# mistyR model training functions
# Copyleft (ɔ) 2020 Jovan Tanevski <jovan.tanevski@uni-heidelberg.de>

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
build_model <- function(views, target, bypass.intra = FALSE, seed = 42,
                        cv.folds = 10, cached = FALSE, ...) {
  
  cache.location <- R.utils::getAbsolutePath(paste0(
    ".misty.temp", .Platform$file.sep,
    views[["misty.uniqueid"]]
  ))

  if (cached && !dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  }

  expr <- views[["intraview"]][["data"]]

  target.vector <- expr %>% dplyr::pull(target)

  # merge ellipsis with default algorithm arguments
  algo.arguments <- list(
    num.trees = 100, importance = "impurity",
    verbose = FALSE, num.threads = 1, seed = seed,
    dependent.variable.name = target
  )

  ellipsis.args <- list(...)
  ellipsis.args.text <- paste(names(ellipsis.args), ellipsis.args,
    sep = ".", collapse = "."
  )

  if (!(length(ellipsis.args) == 0)) {
    algo.arguments <- rlist::list.merge(algo.arguments, ellipsis.args)
  }

  # returns a list of models
  model.views <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map(function(view) {
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

        model.view <- do.call(
          ranger::ranger,
          c(
            list(data = transformed.view.data),
            algo.arguments
          )
        )
        if (cached) {
          readr::write_rds(model.view, model.view.cache.file)
        }
      }

      return(model.view)
    })

  # make oob predictions
  oob.predictions <- model.views %>%
    purrr::map(~ .x$predictions) %>%
    rlist::list.cbind() %>%
    tibble::as_tibble(.name_repair = make.names) %>%
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
    model.views[["intraview"]]$predictions %>%
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

    multi.RMSE <- caret::RMSE(multi.view.prediction, target.vector[test.fold])
    multi.R2 <- caret::R2(multi.view.prediction, target.vector[test.fold],
      formula = "traditional"
    )

    tibble::tibble(
      intra.RMSE = intra.RMSE, intra.R2 = 100*intra.R2,
      multi.RMSE = multi.RMSE, multi.R2 = 100*multi.R2
    )
  })

  final.model <- list(
    meta.model = combined.views,
    model.views = model.views,
    performance.estimate = performance.estimate
  )

  return(final.model)
}
