#' @importFrom magrittr %>% 
#' @importFrom rlang !! :=
.onLoad <- function(libname, pkgname){
  suppressWarnings(future::plan(future::multiprocess))
}


# pass ellipsis to build_model
#' Run MISTy
#'
#' @param views 
#' @param results.folder 
#' @param seed 
#' @param target.subset 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples #TBD
run_misty <- function(views, results.folder = "MVResults",
                                 seed = 42, target.subset = NULL, ...) {
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



  if (!file.exists(paste0(
    results.folder, .Platform$file.sep,
    "coefficients.txt"
  ))) {
    write(header, file = paste0(
      results.folder, .Platform$file.sep,
      "coefficients.txt"
    ))
  } else {
    cat("Coefficients file already exists. Appending!\n")
  }


  header <- "target intra.RMSE intra.R2 multi.RMSE multi.R2 p.RMSE p.R2"

  if (!file.exists(paste0(
    results.folder, .Platform$file.sep,
    "performance.txt"
  ))) {
    write(header, file = paste0(
      results.folder, .Platform$file.sep,
      "performance.txt"
    ))
  } else {
    cat("Performance file already exists. Appending!\n")
  }

  targets <- switch(class(target.subset),
    "numeric" = colnames(expr)[target.subset],
    "integer" = colnames(expr)[target.subset],
    "character" = target.subset,
    "NULL" = colnames(expr),
    NULL
  )

  ranger.available <- require("ranger", quietly = TRUE)

  targets %>% furrr::future_map_chr(function(target) {
    target.model <- build_model(views, target, seed, ...)

    combined.views <- target.model[["meta.model"]]

    model.summary <- summary(combined.views)

    # coefficient values and p-values
    # WARNING: hardcoded column index
    coeff <- c(model.summary$coefficients[, 1], model.summary$coefficients[, 4])

    write(paste(target, paste(coeff, collapse = " ")),
      file = paste0(
        results.folder, .Platform$file.sep,
        "coefficients.txt"
      ), append = TRUE
    )

    # raw importances
    target.model[["model.views"]] %>% purrr::walk2(
      view.abbrev,
      function(model.view, abbrev) {
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
    )

    # performance
    if(sum(target.model[["performance.estimate"]] < 0) > 0){
      warning(paste("Negative performance detected and replaced with 0 for target", target))
    }
    
    performance.estimate <- target.model[["performance.estimate"]] %>%
      mutate_if(~sum(. < 0) > 0, ~pmax(., 0))
    performance.summary <- c(
      performance.estimate %>% colMeans(),
      tryCatch(t.test(performance.estimate %>% dplyr::pull(intra.RMSE),
        performance.estimate %>% dplyr::pull(multi.RMSE),
        alternative = "greater"
      )$p.value, error = function(e){1}),
      tryCatch(t.test(performance.estimate %>% dplyr::pull(intra.R2),
        performance.estimate %>% dplyr::pull(multi.R2),
        alternative = "less"
      )$p.value, error = function(e){1})
    )

    write(paste(target, paste(performance.summary, collapse = " ")),
      file = paste0(results.folder, .Platform$file.sep, "performance.txt"),
      append = TRUE
    )

    return(target)
  }, .progress = TRUE)
}