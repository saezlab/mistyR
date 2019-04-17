# improvement estimation
# pass ellipsis to build_model
estimate_improvement <- function(views, results.folder = "MVResults",
                                 seed = 42, folds = 10, target.subset = NULL, ...) {
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
          predict(model, view[["data"]] %>%
            mutate(!!target := target.vector[fold]), seed = seed)$predictions
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


# pass ellipsis to build_model
estimate_importances <- function(views, results.folder = "MVResults",
                                 seed = 42, target.subset = NULL, ...) {
  if (!dir.exists(results.folder)) dir.create(results.folder, recursive = T)

  view.abbrev <- views %>%
    list.remove(c("misty.uniqueid")) %>%
    map_chr(~ .x[["abbrev"]])


  header <- str_glue("target {views} {p.views}",
    views = paste0(view.abbrev, collapse = " "),
    p.views = paste0("p.", view.abbrev, collapse = " "),
    .sep = " "
  )

  expr <- views[["intracellular"]][["data"]]

  write(header, file = paste0(
    results.folder, .Platform$file.sep,
    "coefficients.txt"
  ))

  targets <- switch(class(target.subset),
    "numeric" = colnames(expr)[target.subset],
    "integer" = colnames(expr)[target.subset],
    "character" = target.subset,
    "NULL" = colnames(expr),
    NULL
  )

  targets %>% future_map_chr(function(target) {
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
    target.model[["model.views"]] %>% walk2(
      view.abbrev,
      function(model.view, abbrev) {
        model.view.imps <- model.view$variable.importance
        imps <- tibble(
          target = names(model.view.imps),
          imp = model.view.imps
        )
        write_csv(imps,
          path = paste0(
            results.folder, .Platform$file.sep,
            "importances_", target, "_", abbrev, ".txt"
          )
        )
      }
    )

    return(target)
  }, .progress = TRUE)
}
