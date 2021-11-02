
library(mistyR)
data("synthetic")
library(tidyverse)

# function to merge arguments
merge_2 <- function(l1, l2) {
  
  n1 <- names(l1)
  n2 <- names(l2)
  diff <- n1[!(n1 %in% n2)]
  n1_list <- diff %>%
    purrr::set_names() %>%
    purrr::map(function(name) l1[[name]])
  
  union <- n2[!(n2 %in% diff)]
  n2_list <- union %>%
    purrr::set_names() %>%
    purrr::map(function(name) l2[[name]])
  return(c(n1_list, n2_list))
}

#####

# keep the ellipsis to pass some parameters to ranger
ranger_model <- function(view_data, target, seed, ...) {
  
  ellipsis.args <- list(...)
  
  # default ranger arguments
  algo.arguments <- list(
    formula = stats::as.formula(paste0(target, " ~ .")),
    data = view_data,
    num.trees = 100,
    importance = "impurity",
    mtry = NULL, 
    verbose = FALSE, 
    num.threads = 1,
    seed = seed)
  
  if (!(length(ellipsis.args) == 0)) {
    algo.arguments <- merge_2(algo.arguments, ellipsis.args)
  }
  
  model <- do.call(ranger::ranger, algo.arguments)
  
  predictions <- tibble::tibble(index = seq_len(nrow(view_data)), 
                                prediction = model$predictions)
  
  list(unbiased.predictions = predictions, 
       importances = model$variable.importance)
}

#####

expr <- synthetic$synthetic1 %>% select(-c(row, col, type))
test <- ranger_model(expr, "ECM", splitrule = "extratrees", 42)
test$unbiased.predictions
test$importances

view_data <- expr

#####

bagged_earth_model <- function(view_data, target, seed, n.vars...) {
  
  ellipsis.args <- list(...)
  
  # how many predictors and how many variables to consider for each bag
  n.vars <- ifelse(is.null(n.vars), ncol(view_data)-1, n.vars)
  predictors <- colnames(view_data)[colnames(view_data) != target]
  
  # generate the bags
  bags <- withr::with_seed(
    seed,
    caret::createResample(1:nrow(view_data), times = 100)
  )
  
  # build one model for each bag, return oob predictions and importances
  models <- purrr::map(bags, function(bag) {
    
    vars <- sample(predictors, n.vars)
    
    algo.arguments <- list(
      formula = stats::as.formula(paste0(target, " ~ .")),
      data = view_data[bag, c(target, vars)]
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_2(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(earth::earth, algo.arguments)
    
    oob <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% bag)]
    
    pred <- predict(model, view_data[oob, vars])
    list(model = model, 
         prediction = tibble::tibble(index = oob, prediction = as.vector(pred)))
  })
  
  # Generate the OOB predictions, using the average for each sample
  predictions <- purrr::map_dfr(models, function(model) {
    tibble::tibble(model$prediction)
  }) %>% 
    dplyr::group_by(index) %>%
    dplyr::summarise(prediction = mean(prediction)) %>%
    dplyr::arrange(index)
  
  assertthat::assert_that(nrow(predictions) == nrow(view_data),
      msg = "There are too few bags to get OOB predictions for all observations.
      Consider increasing the number of bags or using CV.")
  
  importances <- purrr::map_dfr(models, function(model) {
    coefs <- earth::evimp(model$model, trim = FALSE, sqrt. = TRUE)[, 6] 
    names(coefs) <- stringr::str_remove(names(coefs), "-unused")
    coefs
  }) %>% colMeans(na.rm = TRUE)
  
  
  list(unbiased.predictions = predictions, 
       importances = importances)
}

#####

expr <- synthetic$synthetic1 %>% select(-c(row, col, type))
test <- bagged_earth_model(expr, "ECM", 42, NULL)
test$unbiased.predictions
test$importances

#####

bagged_linear_model = function(view_data, target, seed, n.vars...) {
  
  ellipsis.args <- list(...)
  
  # how many predictors and how many variables to consider for each bag
  n.vars <- ifelse(is.null(n.vars), ncol(view_data)-1, n.vars)
  predictors <- colnames(view_data)[colnames(view_data) != target]
  
  # generate the bags
  bags <- withr::with_seed(
    seed,
    caret::createResample(1:nrow(view_data), times = 100)
  )
  
  # build one model for each bag, return oob predictions and importances
  models <- purrr::map(bags, function(bag) {
    
    vars <- sample(predictors, n.vars)
    
    algo.arguments <- list(
      formula = stats::as.formula(paste0(target, " ~ .")),
      data = view_data[bag, c(target, vars)]
    )
    
    if (!(length(ellipsis.args) == 0)) {
      algo.arguments <- merge_2(algo.arguments, ellipsis.args)
    }
    
    model <- do.call(stats::lm, algo.arguments)
    
    oob <- seq.int(1, nrow(view_data))[!(seq.int(1, nrow(view_data)) %in% bag)]
    
    pred <- predict.lm(model, view_data[oob, vars])

    list(model = model, 
         prediction = tibble::tibble(index = oob, prediction = as.vector(pred)))
  })
  
  # Generate the OOB predictions, using the average for each sample
  predictions <- purrr::map_dfr(models, function(model) {
    tibble::tibble(model$prediction)
  }) %>% 
    dplyr::group_by(index) %>%
    dplyr::summarise(prediction = mean(prediction)) %>%
    dplyr::arrange(index)
  
  assertthat::assert_that(nrow(predictions) == nrow(view_data),
      msg = "There are too few bags to get OOB predictions for all observations.
      Consider increasing the number of bags or using CV.")
  
  importances <- purrr::map_dfr(models, function(model) {
    coefs <- purrr::map_dfr(models, function(model) {
      model$model$coefficients }) %>%
      colMeans(na.rm = TRUE)
    coefs[names(coefs) != "(Intercept)"]
  }) %>% colMeans(na.rm = TRUE)
  
  
  list(unbiased.predictions = predictions, 
       importances = importances)
}

#####

expr <- synthetic$synthetic1 %>% select(-c(row, col, type))
test <- bagged_linear_model(expr, "ECM", 42, NULL)
test$unbiased.predictions
test$importances

#####

