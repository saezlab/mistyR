#' @importFrom magrittr %>%
#' @importFrom rlang !! :=
.onLoad <- function(libname, pkgname) {
  suppressWarnings(future::plan(future::multiprocess))
}


# pass ellipsis to build_model

#' Run MISTy
#'
#' @param views
#' @param results.folder
#' @param seed
#' @param target.subset
#' @param cv.folds
#' @param cached
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
run_misty <- function(views, results.folder = "results",
                      seed = 42, target.subset = NULL, cv.folds = 10, cached = TRUE, ...) {
  if (!dir.exists(results.folder)) dir.create(results.folder, recursive = T)

  view.abbrev <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map_chr(~ .x[["abbrev"]])


  header <- stringr::str_glue("target intercept {views} p.intercept {p.views}",
    views = paste0(view.abbrev, collapse = " "),
    p.views = paste0("p.", view.abbrev, collapse = " "),
    .sep = " "
  )

  expr <- views[["intracellular"]][["data"]]



  coef.file <- paste0(results.folder, .Platform$file.sep, "coefficients.txt")
  coef.lock <- paste0(results.folder, .Platform$file.sep, "coefficients.txt.lock")
  on.exit(file.remove(coef.lock))

  if (!file.exists(coef.file)) {
    current.lock <- filelock::lock(coef.lock)
    write(header, file = coef.file)
    filelock::unlock(current.lock)
  } else {
    message("Coefficients file already exists. Appending!\n")
  }


  header <- "target intra.RMSE intra.R2 multi.RMSE multi.R2 p.RMSE p.R2"

  perf.file <- paste0(results.folder, .Platform$file.sep, "performance.txt")
  perf.lock <- paste0(results.folder, .Platform$file.sep, "performance.txt.lock")
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

        readr::write_csv(imps,
          path = paste0(
            results.folder, .Platform$file.sep,
            "importances_", target, "_", abbrev, ".txt"
          )
        )
      }
    )

    # performance
    if (sum(target.model[["performance.estimate"]] < 0) > 0) {
      warning(paste("Negative performance detected and replaced with 0 for target", target))
    }

    performance.estimate <- target.model[["performance.estimate"]] %>%
      dplyr::mutate_if(~ sum(. < 0) > 0, ~ pmax(., 0))
    performance.summary <- c(
      performance.estimate %>% colMeans(),
      tryCatch(t.test(performance.estimate %>% dplyr::pull(intra.RMSE),
        performance.estimate %>% dplyr::pull(multi.RMSE),
        alternative = "greater"
      )$p.value, error = function(e) {
        1
      }),
      tryCatch(t.test(performance.estimate %>% dplyr::pull(intra.R2),
        performance.estimate %>% dplyr::pull(multi.R2),
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
  }, .progress = TRUE, .options = furrr::future_options(seed = TRUE))



  return(results.folder)
}
