build_model <- function(views, target, seed = 42, cached = TRUE) {
  set.seed(seed)
  
  cache.location <- paste0(
    ".misty.temp", .Platform$file.sep,
    views[["misty.uniqueid"]]
  )
  
  if (!dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  }
  
  expr <- views[["intracellular"]][["data"]]
  
  # use the non formula signature of randomForest for memory efficiency
  target.vector <- expr %>% pull(target)
  
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
        target.index <- match(target, colnames(view[["data"]]))
        model.view <- randomForest(
          x = view[["data"]] %>% select(-target.index),
          y = target.vector, ntree = 100
        )
        
        if (cached) {
          write_rds(model.view, model.view.cache.file)
        }
      }
      
      return(model.view)
    })
  
  # make oob predictions
  oob.predictions <- model.views %>%
    map(~ predict(.x)) %>%
    list.cbind() %>%
    as_tibble() %>%
    add_column(!!target := target.vector)
  
  # train lm on above
  combined.views <- lm(as.formula(paste0(target, "~.")), oob.predictions)
  
  # make final.model an object from class misty.model?
  final.model <- list(
    meta.model = combined.views,
    model.views = model.views
  )
  
  return(final.model)
}