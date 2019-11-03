# pass ellipsis to build_model
estimate_importances <- function(views, results.folder = "MVResults",
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

    performance.estimate <- target.model[["performance.estimate"]]
    performance.summary <- c(
      performance.estimate %>% colMeans(),
      t.test(performance.estimate %>% dplyr::pull(intra.RMSE),
        performance.estimate %>% dplyr::pull(multi.RMSE),
        alternative = "greater"
      )$p.value,
      t.test(performance.estimate %>% dplyr::pull(intra.R2),
        performance.estimate %>% dplyr::pull(multi.R2),
        alternative = "less"
      )$p.value
    )

    write(paste(target, paste(performance.summary, collapse = " ")),
      file = paste0(results.folder, .Platform$file.sep, "performance.txt"),
      append = TRUE
    )

    return(target)
  }, .progress = TRUE)
}



# DEPRECATED
# improvement estimation
# pass ellipsis to build_model
estimate_improvement <- function(views, results.folder = "MVResults",
                                 seed = 42, folds = 10, target.subset = NULL, ...) {
  .Deprecated("estimate_importances",
    msg = "the improvements are now estimated at the level of the meta-model"
  )

  if (!dir.exists(results.folder)) dir.create(results.folder, recursive = T)

  expr <- views[["intracellular"]][["data"]]

  header <- "target intra.RMSE intra.R2 multi.RMSE multi.R2"

  write(header, file = paste0(
    results.folder, .Platform$file.sep,
    "performance.txt"
  ))


  targets <- switch(class(target.subset),
    "numeric" = colnames(expr)[target.subset],
    "integer" = colnames(expr)[target.subset],
    "character" = target.subset,
    "NULL" = colnames(expr),
    NULL
  )

  ranger.available <- "ranger" %in% rownames(installed.packages())

  # nested futures! a proper topology should be defined using plan()
  # e.g. plan(list(tweak(multiprocess, workers = 2),
  #               tweak(multiprocess, workers = 4)))
  # for handling 2 targets times 4 folds in parralel
  targets %>% future_map_chr(function(target) {
    set.seed(seed)

    target.vector <- expr %>% pull(target)

    test.folds <- createFolds(seq(nrow(expr)), k = folds)

    performance.metrics <- test.folds %>% future_map_dfr(function(fold) {
      # remove test from views
      train.views <- views %>%
        list.remove(c("misty.uniqueid")) %>%
        map(function(view) {
          list(abbrev = view[["abbrev"]], data = view[["data"]] %>%
            slice(-fold))
        }) %>%
        append(list(misty.uniqueid = views[["misty.uniqueid"]]))

      test.views <- views %>%
        list.remove(c("misty.uniqueid")) %>%
        map(function(view) {
          list(abbrev = view[["abbrev"]], data = view[["data"]] %>% slice(fold))
        })

      model.trained <- build_model(train.views, target, seed, FALSE, ...)


      all.predictions <- model.trained[["model.views"]] %>%
        map2(test.views, function(model, view) {
          if (ranger.available) {
            predict(model, view[["data"]] %>%
              mutate(!!target := target.vector[fold]), seed = seed)$predictions
          } else {
            predict(model, view[["data"]] %>%
              mutate(!!target := target.vector[fold]))
          }
        })


      intra.prediction <- all.predictions$intracellular
      multi.view.prediction <- predict(
        model.trained[["meta.model"]],
        as_tibble(all.predictions) %>% mutate(!!target := target.vector[fold])
      )


      intra.RMSE <- RMSE(intra.prediction, target.vector[fold])
      intra.R2 <- R2(intra.prediction, target.vector[fold],
        formula = "traditional"
      )

      multi.RMSE <- RMSE(multi.view.prediction, target.vector[fold])
      multi.R2 <- R2(multi.view.prediction, target.vector[fold],
        formula = "traditional"
      )

      tibble(
        intra.RMSE = intra.RMSE, intra.R2 = intra.R2,
        multi.RMSE = multi.RMSE, multi.R2 = multi.R2
      )
    })

    performance.metrics.summary <- performance.metrics %>% colMeans()

    write(paste(target, paste(performance.metrics.summary, collapse = " ")),
      file = paste0(results.folder, .Platform$file.sep, "performance.txt"),
      append = TRUE
    )

    return(target)
  }, .progress = TRUE)
}
