#' @importFrom magrittr %>%
#' @importFrom rlang !! :=
.onLoad <- function(libname, pkgname) {
  suppressWarnings(future::plan(future::multiprocess))
}

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
  ranger.available <- require("ranger", quietly = TRUE)
  if (!dir.exists(results.folder)) dir.create(results.folder, recursive = T)

  view.abbrev <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map_chr(~ .x[["abbrev"]])

  coefficient.file <- initiate_coefficient_file(view.abbrev, results.folder)
  performance.file <- initiate_performance_file(results.folder)
  targets <- get_targets_names(views, target.subset)
  
  message("Training models")
  targets %>% furrr::future_map_chr(function(target, ...) {
    target.model <- build_model(views, target, seed, cv.folds, cached, ...)
    
    combined.views <- target.model[["meta.model"]]
    model.summary <- summary(combined.views)
  
    # coefficient values and p-values
    # WARNING: hardcoded column index
    coeff <- c(model.summary$coefficients[, 1], model.summary$coefficients[, 4])
    coefficient.data <- paste(target, paste(coeff, collapse = " "))
    write_results_to_file(coefficient.file, coefficient.data)
    
    target.model[["model.views"]] %>% purrr::walk2(
      .y=view.abbrev,
      .f=~calculate_raw_importances(., .y, target, results.folder, ranger.available)
    )

    performance.summary <- estimate_performance(target.model[["performance.estimate"]], target)
    performance.data <- paste(target, paste(performance.summary, collapse = " "))
    write_results_to_file(performance.file, performance.data)

    return(target)
  }, .progress = TRUE)

  return(results.folder)
}


initiate_coefficient_file <- function(view.abbrev, results.folder) {
  
  header <- stringr::str_glue("target intercept {views} p.intercept {p.views}",
                              views = paste0(view.abbrev, collapse = " "),
                              p.views = paste0("p.", view.abbrev, collapse = " "),
                              .sep = " ")
  
  coef.file <- paste0(results.folder, .Platform$file.sep, "coefficients.txt")
  
  if (!file.exists(coef.file)) {
    write_results_to_file(coef.file, header)
  } else {
    message("Coefficients file already exists. Appending!\n")
  }
  
  return(coef.file)
}


initiate_performance_file <- function(results.folder) {
  
  header <- "target intra.RMSE intra.R2 multi.RMSE multi.R2 p.RMSE p.R2"
  
  perf.file <- paste0(results.folder, .Platform$file.sep, "performance.txt")
  
  if (!file.exists(perf.file)) {
    write_results_to_file(perf.file, header)
  } else {
    message("Performance file already exists. Appending!\n")
  }
  
  return(perf.file)
}

write_results_to_file <- function(file, data) {
  current.lock <- filelock::lock(file)
  write(data, file = file, append = TRUE)
  filelock::unlock(current.lock)

}

get_targets_names <- function(views, target.subset = NULL) {
  all.target.names <- colnames( views[["intracellular"]][["data"]] )
  
  target.names <- switch(class(target.subset),
         "numeric" = all.target.names[target.subset],
         "integer" = all.target.names[target.subset],
         "character" = target.subset,
         "NULL" = all.target.names,
         NULL
  )
  
  return(target.names)
}

calculate_raw_importances <- function(model.view, abbrev, target, results.folder, ranger.available=FALSE) {
  if (ranger.available) {
    model.view.imps <- model.view$variable.importance
    targets <- names(model.view.imps)
  } else {
    model.view.imps <- randomForest::importance(model.view, type = 2)
    targets <- rownames(model.view.imps)
  }
  
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

estimate_performance <- function(performance.estimate, target) {
  if (sum(performance.estimate < 0) > 0) {
    warning(paste("Negative performance detected and replaced with 0 for target", target))
  }
  
  performance.estimate <- performance.estimate %>%
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
  return(performance.summary)
}