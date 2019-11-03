# elipsis passed to ranger or randomForest
build_model <- function(views, target, seed = 42, cached = TRUE, cv.folds = 10, ...) {
  set.seed(seed)

  cache.location <- paste0(
    ".misty.temp", .Platform$file.sep,
    views[["misty.uniqueid"]]
  )

  if (!dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  }

  expr <- views[["intracellular"]][["data"]]

  target.vector <- expr %>% dplyr::pull(target)

  ranger.available <- require("ranger", quietly = TRUE)

  # merge ellipsis with default algorithm arguments
  if (ranger.available) {
    algo.arguments <- list(
      num.trees = 100, importance = "impurity",
      verbose = FALSE, num.threads = 1, seed = seed,
      dependent.variable.name = target
    )
  } else {
    algo.arguments <- list(ntree = 100)
  }

  if (!length(list(...)) == 0) {
    algo.arguments <- rlist::list.merge(algo.arguments, list(...))
  }

  # returns a list of models
  model.views <- views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map(function(view) {
      model.view.cache.file <-
        paste0(
          cache.location, .Platform$file.sep,
          "model.", view[["abbrev"]], ".", target, ".rds"
        )

      if (file.exists(model.view.cache.file) & cached) {
        model.view <- readr::read_rds(model.view.cache.file)
      } else {
        if (ranger.available) {
          model.view <- do.call(
            ranger::ranger,
            c(
              list(data = (view[["data"]] %>%
                dplyr::mutate(!!target := target.vector))),
              algo.arguments
            )
          )
        } else {
          target.index <- match(target, colnames(view[["data"]]))
          model.view <- do.call(
            randomForest::randomForest,
            c(
              list(
                x = view[["data"]] %>% dplyr::select(-target.index),
                y = target.vector
              ), algo.arguments
            )
          )
        }
        if (cached) {
          readr::write_rds(model.view, model.view.cache.file)
        }
      }

      return(model.view)
    })

  # make oob predictions
  oob.predictions <- model.views %>%
    purrr::map(~ if (ranger.available) {
      .x$predictions
    } else {
      .x$predicted
    }) %>%
    rlist::list.cbind() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(!!target := target.vector)

  # train lm on above
  combined.views <- lm(as.formula(paste0(target, "~.")), oob.predictions)

  # cv performance estimate
  test.folds <- caret::createFolds(target.vector, k = cv.folds)

  intra.view.only <- (if (ranger.available) {
    model.views[["intracellular"]]$predictions
  }
  else {
    model.views[["intracellular"]]$predicted
  }) %>%
    tibble::enframe(name = NULL) %>%
    dplyr::mutate(!!target := target.vector)

  performance.estimate <- test.folds %>% purrr::map_dfr(function(test.fold) {
    meta.intra <- lm(
      as.formula(paste0(target, "~.")),
      intra.view.only %>% dplyr::slice(-test.fold)
    )
    meta.multi <- lm(
      as.formula(paste0(target, "~.")),
      oob.predictions %>% dplyr::slice(-test.fold)
    )

    intra.prediction <- predict(meta.intra, intra.view.only %>%
      dplyr::slice(test.fold))
    multi.view.prediction <- predict(meta.multi, oob.predictions %>%
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
      intra.RMSE = intra.RMSE, intra.R2 = intra.R2,
      multi.RMSE = multi.RMSE, multi.R2 = multi.R2
    )
  })


  # make final.model an object from class misty.model?
  final.model <- list(
    meta.model = combined.views,
    model.views = model.views,
    performance.estimate = performance.estimate
  )

  return(final.model)
}
