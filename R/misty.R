# MISTy runner
# Copyright (c) 2020 Jovan Tanevski [jovan.tanevski@uni-heidelberg.de]


#' @importFrom magrittr %>%
#' @importFrom rlang !! := .data
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("mistyR is able to run computationally intensive functions
  in parallel. Please consider specifying a future::plan(). For example by running
  future::plan(future::multisession) before calling mistyR functions.")
}

#' Train MISTy models
#'
#' Trains multi-view models for all target markers, estimates the performance,
#' the contributions of the view specific models and the importance of predictor
#' markers for each target marker.
#'
#' Default values passed to \code{\link[ranger]{ranger}()} for training the
#' view-specific models: \code{num.trees = 100}, \code{importance = "impurity"},
#' \code{num.threads = 1}, \code{seed = seed}.
#'
#' @param views view composition.
#' @param results.folder path to the top level folder to store raw results.
#' @param seed seed used for random sampling to ensure reproducibility.
#' @param target.subset subset of targets to train models for. If \code{NULL},
#'     models will be trained for markers in the intraview.
#' @param cv.folds number of cross-validation folds to consider for estimating
#'     the performance of the multi-view models.
#' @param cached a \code{logical} indicating whether to cache the trained models
#'     and to reuse previously cached ones if they already exist for this sample.
#' @param ... all additional parameters are passed to
#'     \code{\link[ranger]{ranger}()} for training the view-specific models
#'     (see Details for defaults).
#'
#' @return Path to the results folder that can be passed to
#'     \code{\link{collect_results}()}.
#'
#' @seealso \code{\link{create_initial_view}()} for
#'     starting a view composition.
#'
#' @examples
#' # Create a view composition of an intraview and a paraview with radius 10 then
#' # run MISTy for a single sample.
#'
#' library(dplyr)
#'
#' # get the expression data
#' data("synthetic")
#' expr <- synthetic[[1]] %>% select(-c(row, col, type))
#' # get the coordinates for each cell
#' pos <- synthetic[[1]] %>% select(row, col)
#'
#' # compose
#' misty.views <- create_initial_view(expr) %>% add_paraview(pos, l = 10)
#'
#' # run with default parameters
#' run_misty(misty.views)
#'
#' # Alternatives
#' \dontrun{
#'
#' create_initial_view(expr) %>%
#'   add_paraview(pos, l = 10) %>%
#'   run_misty()
#' }
#' @export
run_misty <- function(views, results.folder = "results", seed = 42,
                      target.subset = NULL, cv.folds = 10, cached = TRUE, ...) {
  normalized.results.folder <- normalizePath(results.folder)

  if (!dir.exists(normalized.results.folder)) {
    dir.create(normalized.results.folder, recursive = TRUE)
  }

  view.abbrev <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map_chr(~ .x[["abbrev"]])


  header <- stringr::str_glue("target intercept {views} p.intercept {p.views}",
    views = paste0(view.abbrev, collapse = " "),
    p.views = paste0("p.", view.abbrev, collapse = " "),
    .sep = " "
  )

  expr <- views[["intraview"]][["data"]]

  assertthat::assert_that(nrow(expr) >= cv.folds,
    msg = "The data has less rows than the requested number of cv folds."
  )


  coef.file <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "coefficients.txt"
  )
  coef.lock <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "coefficients.txt.lock"
  )
  on.exit(file.remove(coef.lock))

  if (!file.exists(coef.file)) {
    current.lock <- filelock::lock(coef.lock)
    write(header, file = coef.file)
    filelock::unlock(current.lock)
  } else {
    message("Coefficients file already exists. Appending!\n")
  }


  header <- "target intra.RMSE intra.R2 multi.RMSE multi.R2 p.RMSE p.R2"

  perf.file <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "performance.txt"
  )
  perf.lock <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "performance.txt.lock"
  )
  on.exit(file.remove(perf.lock), add = TRUE)

  if (!file.exists(perf.file)) {
    current.lock <- filelock::lock(perf.lock)
    write(header, file = perf.file)
    filelock::unlock(current.lock)
  } else {
    message("Performance file already exists. Appending!\n")
  }

  targets <- switch(class(target.subset),
    "numeric" = colnames(expr)[target.subset],
    "integer" = colnames(expr)[target.subset],
    "character" = target.subset,
    "NULL" = colnames(expr),
    NULL
  )

  message("Training models")
  targets %>% furrr::future_map_chr(function(target, ...) {
    target.model <- build_model(views, target, seed, cv.folds, cached, ...)

    combined.views <- target.model[["meta.model"]]

    model.summary <- summary(combined.views)

    # coefficient values and p-values
    # WARNING: hardcoded column index
    coeff <- c(model.summary$coefficients[, 1], model.summary$coefficients[, 4])

    current.lock <- filelock::lock(coef.lock)
    write(paste(target, paste(coeff, collapse = " ")),
      file = coef.file, append = TRUE
    )
    filelock::unlock(current.lock)

    # raw importances
    target.model[["model.views"]] %>% purrr::walk2(
      view.abbrev,
      function(model.view, abbrev) {
        model.view.imps <- model.view$variable.importance
        targets <- names(model.view.imps)

        imps <- tibble::tibble(
          target = targets,
          imp = model.view.imps
        )

        readr::write_csv(
          imps,
          paste0(
            normalized.results.folder, .Platform$file.sep,
            "importances_", target, "_", abbrev, ".txt"
          )
        )
      }
    )

    # performance
    if (sum(target.model[["performance.estimate"]] < 0) > 0) {
      warning.message <-
        paste(
          "Negative performance detected and replaced with 0 for target",
          target
        )
      warning(warning.message)
    }

    performance.estimate <- target.model[["performance.estimate"]] %>%
      dplyr::mutate_if(~ sum(. < 0) > 0, ~ pmax(., 0))
    performance.summary <- c(
      performance.estimate %>% colMeans(),
      tryCatch(stats::t.test(performance.estimate %>% dplyr::pull(.data$intra.RMSE),
        performance.estimate %>% dplyr::pull(.data$multi.RMSE),
        alternative = "greater"
      )$p.value, error = function(e) {
        1
      }),
      tryCatch(stats::t.test(performance.estimate %>% dplyr::pull(.data$intra.R2),
        performance.estimate %>% dplyr::pull(.data$multi.R2),
        alternative = "less"
      )$p.value, error = function(e) {
        1
      })
    )

    current.lock <- filelock::lock(perf.lock)
    write(paste(target, paste(performance.summary, collapse = " ")),
      file = perf.file, append = TRUE
    )
    filelock::unlock(current.lock)

    return(target)
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))



  return(normalized.results.folder)
}
