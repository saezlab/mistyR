# MISTy runner
# Copyleft (ɔ) 2020 Jovan Tanevski [jovan.tanevski@uni-heidelberg.de]

#' @importFrom rlang !! := .data
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("mistyR is able to run computationally intensive functions
  in parallel. Please consider specifying a future::plan(). For example by running
  future::plan(future::multisession) before calling mistyR functions.")
}

#' @importFrom dplyr %>%
#' @inherit run_misty examples
#' @export
dplyr::`%>%`

#' Train MISTy models
#'
#' Trains multi-view models for all target markers, estimates the performance,
#' the contributions of the view specific models and the importance of predictor
#' markers for each target marker.
#'
#' If \code{bypass.intra} is set to \code{TRUE} all variable in the intraview
#' the intraview data will be treated as targets only. The baseline intraview
#' model in this case is a trivial model that predicts the average of each
#' target. If the intraview has only one variable this switch is automatically
#' set to \code{TRUE}.
#'
#' Default ML model to train the view-specific views is 
#' \code{\link[ranger]{ranger}()}. The following parameters are the default
#' configuration: \code{num.trees = 100}, \code{importance = "impurity"},
#' \code{num.threads = 1}, \code{seed = seed}.
#'
#' @param views view composition.
#' @param results.folder path to the top level folder to store raw results.
#' @param seed seed used for random sampling to ensure reproducibility.
#' @param target.subset subset of targets to train models for. If \code{NULL},
#'     models will be trained for markers in the intraview.
#' @param bypass.intra a \code{logical} indicating whether to train a baseline
#'     model using the intraview data (see Details).
#' @param cv.folds number of cross-validation folds to consider for estimating
#'     the performance of the multi-view models and the view-specific models
#'     if \code{method = "cv"}
#' @param cached a \code{logical} indicating whether to cache the trained models
#'     and to reuse previously cached ones if they already exist for this sample.
#' @param append a \code{logical} indicating whether to append the performance
#'     and coefficient files in the \code{results.folder}. Consider setting to
#'     \code{TRUE} when rerunning a workflow with different \code{target.subset}
#'     parameters.
#' @param method a string indicating whether to use bagging "bag" or cross-
#' valdiation to train the view-specific models and to get unbiased predictions
#' for the linear meta model
#' @param learner a string indicating which ML model to use to model the 
#' views, possible values are "ranger" for random forest, "lm" for linear
#' model, "linearSVM" for support vector machine with linear kernel, and
#' "earth" for multivariate adaptive regression splines.
#' @param n.vars number indicating how many variables should be used for 
#' training the view-specific models if \code{method = "bag"}. For each base
#' learner a sample of size \code{n.vars} from the predictors is taken. 
#' If \code{learner = "ranger"} \code{n.vars} corresponds to \code{mtry}, 
#' meaning the number of variables that are considered at each split.
#' @param n.learners number indicating how many learners should be used for 
#' the bagging model, only relevant if \code{method = "bag"}
#' @param ... all additional parameters are passed to the chosen ML model for
#' training the view-specific models 
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
#' @export
run_misty <- function(views, results.folder = "results", seed = 42,
                      target.subset = NULL, bypass.intra = FALSE, cv.folds = 10,
                      cached = FALSE, append = FALSE, method = "bag",
                      learner = "ranger", n.vars = NULL, n.learners = 100, ...) {
  
  assertthat::assert_that(method %in% c("bag", "cv"),
    msg = "The selected method has to be 'bag' (bagging) or 'cv'
    (cross validation)")
  
  supported.models <- c("ranger", "lm", "svmLinear", "earth")
  assertthat::assert_that(learner %in% supported.models,
    msg = paste0("The selected learner (model) is not supported. Currently, the 
                 following models are supported: ", toString(supported.models)))
  
  normalized.results.folder <- R.utils::getAbsolutePath(results.folder)

  if (!dir.exists(normalized.results.folder)) {
    dir.create(normalized.results.folder, recursive = TRUE)
  }

  on.exit(sweep_cache())

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
  
  if (!(is.null(n.vars))) {
    assertthat::assert_that(n.vars <= (ncol(expr)-1),
      msg = "The number of variables selected must be smaller or equal than the 
      number of column of the input minus 1 or NULL"
    )
  }

  if (ncol(expr) == 1) bypass.intra <- TRUE

  target.var <- apply(expr, 2, stats::sd, na.rm = TRUE)

  assertthat::assert_that(!any(target.var == 0),
    msg = paste(
      "Targets",
      paste(names(which(target.var == 0)),
            collapse = ", "
      ),
      "have zero variance."
    ) )
  
  target.unique <- colnames(expr) %>%
    purrr::set_names() %>%
    purrr::map_int(~ length(unique(expr %>% dplyr::pull(.x))))

  # shouldn't we only check for this requirement if we do not bypass intra?
  # (otherwise if we do one-hot encoding of celltypes we can only use
  # 2 CV folds which may be suboptimal)
  assertthat::assert_that(all(target.unique >= cv.folds),
      msg = paste(
        "Targets",
        paste(names(which(target.unique < cv.folds)),
              collapse = ", "
        ),
        "have fewer unique values than cv.folds"
      )
  )

  coef.file <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "coefficients.txt"
  )
  coef.lock <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "coefficients.txt.lock"
  )
  on.exit(file.remove(coef.lock), add = TRUE)

  if (!append) {
    current.lock <- filelock::lock(coef.lock)
    write(header, file = coef.file)
    filelock::unlock(current.lock)
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

  if (!append) {
    current.lock <- filelock::lock(perf.lock)
    write(header, file = perf.file)
    filelock::unlock(current.lock)
  }

  targets <- switch(class(target.subset),
    "numeric" = colnames(expr)[target.subset],
    "integer" = colnames(expr)[target.subset],
    "character" = target.subset,
    "NULL" = colnames(expr),
    NULL
  )

  message("\nTraining models")
  targets %>% furrr::future_map_chr(function(target, ...) {
    target.model <- build_model(views = views, target = target, method = method, 
                                learner = learner, n.vars = n.vars, 
                                n.learners = n.learners, cv.folds = cv.folds, 
                                bypass.intra = bypass.intra, seed = seed, 
                                cached = cached, ...)

    combined.views <- target.model[["meta.model"]]

    model.summary <- summary(combined.views)

    # coefficient values and p-values
    # WARNING: hardcoded column index
    coeff <- c(
      if (bypass.intra) 0, stats::coef(combined.views),
      if (bypass.intra) 1, model.summary$coefficients[, 4]
    )

    current.lock <- filelock::lock(coef.lock)
    write(paste(target, paste(coeff, collapse = " ")),
      file = coef.file, append = TRUE
    )
    filelock::unlock(current.lock)

    # raw importances
    target.model[["model.importances"]] %>% purrr::walk2(
      view.abbrev,
      function(model.importance, abbrev) {
        targets <- names(model.importance)

        imps <- tibble::tibble(
          target = targets,
          imp = model.importance
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
      tryCatch(stats::t.test(performance.estimate %>%
        dplyr::pull(.data$intra.RMSE),
      performance.estimate %>%
        dplyr::pull(.data$multi.RMSE),
      alternative = "greater"
      )$p.value, error = function(e) {
        warning.message <- paste(
          "t-test of RMSE performance failed with error:",
          e$message
        )
        warning(warning.message)
        1
      }),
      tryCatch(stats::t.test(performance.estimate %>%
        dplyr::pull(.data$intra.R2),
      performance.estimate %>%
        dplyr::pull(.data$multi.R2),
      alternative = "less"
      )$p.value, error = function(e) {
        warning.message <- paste(
          "t-test of R2 performance failed with error:",
          e$message
        )
        warning(warning.message)
        1
      })
    )

    current.lock <- filelock::lock(perf.lock)
    write(paste(target, paste(performance.summary, collapse = " ")),
      file = perf.file, append = TRUE
    )
    filelock::unlock(current.lock)

    return(target)
  }, ..., .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))

  return(normalized.results.folder)
}
