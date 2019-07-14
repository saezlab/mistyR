# elipsis passed to ranger or randomForest
build_model <- function(views, target, seed = 42, cached = TRUE, ...) {
  set.seed(seed)

  cache.location <- paste0(
    ".misty.temp", .Platform$file.sep,
    views[["misty.uniqueid"]]
  )

  if (!dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  }

  expr <- views[["intracellular"]][["data"]]

  target.vector <- expr %>% pull(target)

  ranger.available <- "ranger" %in% rownames(installed.packages())

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
  
  if(!length(list(...))==0)
    algo.arguments <- list.merge(algo.arguments, list(...))

  # returns a list of models
  model.views <- views %>%
    list.remove(c("misty.uniqueid")) %>%
    map(function(view) {
      model.view.cache.file <-
        paste0(
          cache.location, .Platform$file.sep,
          "model.", view[["abbrev"]], ".", target, ".rds"
        )

      if (file.exists(model.view.cache.file) & cached) {
        model.view <- read_rds(model.view.cache.file)
      } else {
        if (ranger.available) {
          model.view <- do.call(
            ranger,
            c(
              list(data = (view[["data"]] %>%
                mutate(!!target := target.vector))),
              algo.arguments
            )
          )
        } else {
          target.index <- match(target, colnames(view[["data"]]))
          model.view <- do.call(
            randomForest,
            c(
              list(
                x = view[["data"]] %>% select(-target.index),
                y = target.vector
              ), algo.arguments
            )
          )
        }
        if (cached) {
          write_rds(model.view, model.view.cache.file)
        }
      }

      return(model.view)
    })

  # make oob predictions
  oob.predictions <- model.views %>%
    map(~if(ranger.available){ .x$predictions } else { predict(.x) }) %>%
    list.cbind() %>%
    as_tibble() %>%
    mutate(!!target := target.vector)
  # train lm on above
  combined.views <- lm(as.formula(paste0(target, "~.-1")), oob.predictions)

  # make final.model an object from class misty.model?
  final.model <- list(
    meta.model = combined.views,
    model.views = model.views
  )

  return(final.model)
}
