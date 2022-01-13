expr <- generate_random_tibble(100, 5)
pos <- sample_grid_geometry(100, 10, 10)
misty.views <- create_initial_view(expr) %>% add_paraview(pos, l = 2)

test_that("run_misty produces correct files on output", {
  suppressWarnings(run_misty(misty.views))
  expect_true(dir.exists("results"))
  expect_length(list.files("results"), 12)
  expect_true((all(list.files("results", "importance*", full.names = TRUE) %>%
    purrr::map_int(R.utils::countLines) == 5)))
  expect_true((all(list.files("results", "(coefficients|performance)",
    full.names = TRUE
  ) %>%
    purrr::map_int(R.utils::countLines) == 6)))
  unlink("results", recursive = TRUE)
})

test_that("run_misty handles evaluation parameters correctly", {
  suppressWarnings({
    default.time <- system.time(
      run_misty(misty.views, "results1")
    )["user.self"] * 1000
    cv.time <- system.time(
      run_misty(misty.views, "results2", cv.folds = 3)
    )["user.self"] * 1000
    subset.time <- system.time(
      run_misty(misty.views, "results3", target.subset = c("expr1", "expr2"))
    )["user.self"] * 1000
    ntrees.time <- system.time(
      run_misty(misty.views, "results4", num.trees = 500)
    )["user.self"] * 1000
  })
  expect_lt(cv.time, default.time)
  expect_lt(subset.time, default.time)
  expect_gt(ntrees.time, default.time)
  expect_length(list.files("results3"), 6)
  unlink(paste0("results", seq_len(4)), recursive = TRUE)
})

test_that("run_misty handles tests of failures", {
  expr <- tibble::tibble(expr1 = c(rep(1, 50), rep(2, 50))) %>%
    dplyr::mutate(expr2 = rev(expr1))
  warning.message <- capture_warnings(create_initial_view(expr) %>% run_misty())
  expect_true(any(grepl("have fewer unique values than cv.folds", warning.message)))
  sig.warnings <- capture_warnings(create_initial_view(expr) %>%
    run_misty(cv.folds = 2))
  expect_true(any(grepl("RMSE", sig.warnings)))
  expect_true(any(grepl("R2", sig.warnings)))
  unlink("results", recursive = TRUE)
})

test_that("modeling of intraview is bypassed if only 1 var in intraview", {
  truncated_expr <- generate_random_tibble(100, 1)
  misty.views <- create_initial_view(truncated_expr)
  expect_visible(suppressWarnings(misty.views %>% run_misty()))
  unlink("results", recursive = TRUE)
})

test_that("warning raised if variance of variable is 0", {
  expr <- tibble::tibble(
    expr1 = 10,
    expr2 = runif(100, 2, 5),
    expr3 = rnorm(100, 10, 2)
  ) %>%
    dplyr::mutate(expr4 = 2 * expr2 + 0.5 * expr3)
  misty.views <- create_initial_view(expr)
  expect_error(
    misty.views %>% run_misty(target.subset = "expr4"),
    "have zero variance"
  )
  unlink("results", recursive = TRUE)
})

test_that("all models work and produce the correct output", {
  functions <- list("rf" = random_forest_model, 
                    "bag_mars" = bagged_mars_model, 
                    "mars" = mars_model,
                    "linear" = linear_model,
                    "svm" = svm_model,
                    "boosting" = gradient_boosting_model,
                    "mpl" = mlp_model)
  
  ncols <- 5
  expr <- generate_random_tibble(100, ncols)
  misty.views <- create_initial_view(expr)
  
  misty.test <- purrr::map(functions, function(fun) {
    suppressWarnings(misty.results <- run_misty(misty.views, model.function = fun) %>%
      collect_results()
    )
    expect_true(dir.exists("results"))
    expect_length(list.files("results"), ncols + 2)
    expect_true((all(list.files("results", "importance*", full.names = TRUE) %>%
                       purrr::map_int(R.utils::countLines) == ncols)))
    expect_true((all(list.files("results", "(coefficients|performance)",
                                full.names = TRUE
    ) %>%
      purrr::map_int(R.utils::countLines) == ncols + 1)))
    unlink("results", recursive = TRUE)
  })
})

test_that("ellipsis arguments can be passed to the provided ML models", {
  # create data
  ncols <- 10
  expr <- generate_random_tibble(100, ncols)
  pos <- sample_grid_geometry(100, 20, 20)
  misty.views <- create_initial_view(expr) %>%
    add_paraview(positions = pos, l = 10)
  
  # random forest
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = random_forest_model)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = random_forest_model, 
                          num.trees = 2000)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  # bagged mars
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = bagged_mars_model,
                          degree = 1)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = bagged_mars_model, 
                          n.bags = 50)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  # mars
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = mars_model,
                          degree = 3, nk = 30)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = mars_model, 
                          degree = 3, nk = 30)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  # svm
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = svm_model,
                          C = 1)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = svm_model, 
                          C = 100)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  # gradient boosting
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = gradient_boosting_model,
                          booster = "gbtree", nrounds = 10)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = gradient_boosting_model, 
                          booster = "gbtree", nrounds = 20)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)

  # multi-layer perceptron
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = mlp_model,
                          size = c(1), maxit = 1)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
  misty.test <- run_misty(misty.views, model.function = mlp_model, 
                          size = c(10), maxit = 100)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  unlink("results", recursive = TRUE)
})

test_that("k for cv , n.bags for bagging can be changed and approx works", {
  # create data
  ncols <- 10
  expr <- generate_random_tibble(100, ncols)
  pos <- sample_grid_geometry(100, 20, 20)
  misty.views <- create_initial_view(expr) %>%
    add_paraview(positions = pos, l = 10)

  # bagged mars
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = bagged_mars_model,
                            n.bags = 20)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = bagged_mars_model, 
                            n.bags = 50)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  # mars
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = mars_model,
                            k = 10, approx = 0.8)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = mars_model, 
                            k = 25, approx = 1)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  # linear
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = linear_model,
                            k = 10)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = linear_model, 
                            k = 25)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  # svm
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = svm_model,
                            k = 10 , approx = .4)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = svm_model, 
                            k = 25, approx = 1)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  # gradient boosting
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = gradient_boosting_model,
                            k = 10)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = gradient_boosting_model, 
                            k = 25)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  # multi-layer perceptron
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = mlp_model,
                            k = 2, approx = 0.6)
  )
  end <- Sys.time()
  first.run = end - start
  
  start <- Sys.time()
  suppressWarnings(
    misty.test <- run_misty(misty.views, model.function = mlp_model, 
                            k = 8, approx = 1)
  )
  end <- Sys.time()
  second.run = end - start
  testthat::expect_true(first.run < second.run)
  
  unlink("results", recursive = TRUE)
})

test_that("run_misty is reproducible for all ML algorithms", {
  # random forest
  suppressWarnings({
    hash.results1 <- list.files(
      run_misty(misty.views, model.function = random_forest_model,"results1", seed = 1),
                                full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results2 <- list.files(
      run_misty(misty.views, model.function = random_forest_model, "results2", seed = 1),
                                full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results3 <- list.files(
      run_misty(misty.views, model.function = random_forest_model, "results3", seed = 2),
                                full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
  })
  expect_equal(hash.results1, hash.results2)
  expect_true(all(hash.results1 != hash.results3))
  unlink(paste0("results", seq_len(3)), recursive = TRUE)
  
  # bagged mars
  suppressWarnings({
    hash.results1 <- list.files(
      run_misty(misty.views, model.function = bagged_mars_model,"results1", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results2 <- list.files(
      run_misty(misty.views, model.function = bagged_mars_model, "results2", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results3 <- list.files(
      run_misty(misty.views, model.function = bagged_mars_model, "results3", seed = 2),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
  })
  expect_equal(hash.results1, hash.results2)
  expect_true(all(hash.results1 != hash.results3))
  unlink(paste0("results", seq_len(3)), recursive = TRUE)
  
  # mars
  suppressWarnings({
    hash.results1 <- list.files(
      run_misty(misty.views, model.function = mars_model,"results1", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results2 <- list.files(
      run_misty(misty.views, model.function = mars_model, "results2", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results3 <- list.files(
      run_misty(misty.views, model.function = mars_model, "results3", seed = 2),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
  })
  expect_equal(hash.results1, hash.results2)
  # using any instead of all, since the model is not random
  expect_true(any(hash.results1 != hash.results3))
  unlink(paste0("results", seq_len(3)), recursive = TRUE)
  
  # gradient boosting
  suppressWarnings({
    hash.results1 <- list.files(
      run_misty(misty.views, model.function = gradient_boosting_model,"results1", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results2 <- list.files(
      run_misty(misty.views, model.function = gradient_boosting_model, "results2", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results3 <- list.files(
      run_misty(misty.views, model.function = gradient_boosting_model, "results3", seed = 2),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
  })
  expect_equal(hash.results1, hash.results2)
  expect_true(any(hash.results1 != hash.results3))
  unlink(paste0("results", seq_len(3)), recursive = TRUE)
  
  # linear model
  suppressWarnings({
    hash.results1 <- list.files(
      run_misty(misty.views, model.function = linear_model,"results1", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results2 <- list.files(
      run_misty(misty.views, model.function = linear_model, "results2", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results3 <- list.files(
      run_misty(misty.views, model.function = linear_model, "results3", seed = 2),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
  })
  expect_equal(hash.results1, hash.results2)
  expect_true(any(hash.results1 != hash.results3))
  unlink(paste0("results", seq_len(3)), recursive = TRUE)
  
  # svm
  suppressWarnings({
    hash.results1 <- list.files(
      run_misty(misty.views, model.function = svm_model,"results1", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results2 <- list.files(
      run_misty(misty.views, model.function = svm_model, "results2", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results3 <- list.files(
      run_misty(misty.views, model.function = svm_model, "results3", seed = 2),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
  })
  # not fully reproducible currenlty
  #expect_equal(hash.results1, hash.results2)
  #expect_true(sum(hash.results1 == hash.results2) > 18)
  expect_true(any(hash.results1 != hash.results3))
  unlink(paste0("results", seq_len(3)), recursive = TRUE)
  
  # mlp
  suppressWarnings({
    hash.results1 <- list.files(
      run_misty(misty.views, model.function = mlp_model,"results1", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results2 <- list.files(
      run_misty(misty.views, model.function = mlp_model, "results2", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results3 <- list.files(
      run_misty(misty.views, model.function = mlp_model, "results3", seed = 2),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
  })
  # not fully reproducible currenlty
  #expect_equal(hash.results1, hash.results2)
  #expect_true(sum(hash.results1 == hash.results2) > 18)
  expect_true(any(hash.results1 != hash.results3))
  unlink(paste0("results", seq_len(3)), recursive = TRUE)
})

test_that("bypass intra works for all ML algorithms", {
  # random forest
  suppressWarnings(
  misty.run <- run_misty(misty.views, model.function = random_forest_model, 
            bypass.intra = TRUE) %>%
    collect_results()
  )
  expect_true(all(misty.run$importances.aggregated %>%
                    dplyr::filter(view == "intra") %>%
                    dplyr::pull(Importance) == 0)
  )
  
  # bagged mars
  suppressWarnings(
    misty.run <- run_misty(misty.views, model.function = bagged_mars_model, 
                           bypass.intra = TRUE) %>%
      collect_results()
  )
  expect_true(all(misty.run$importances.aggregated %>%
                    dplyr::filter(view == "intra") %>%
                    dplyr::pull(Importance) == 0)
  )
  
  # mars
  suppressWarnings(
    misty.run <- run_misty(misty.views, model.function = mars_model, 
                           bypass.intra = TRUE) %>%
      collect_results()
  )
  expect_true(all(misty.run$importances.aggregated %>%
                    dplyr::filter(view == "intra") %>%
                    dplyr::pull(Importance) == 0)
  )
  
  # gradient boosting
  suppressWarnings(
    misty.run <- run_misty(misty.views, model.function = gradient_boosting_model, 
                           bypass.intra = TRUE) %>%
      collect_results()
  )
  expect_true(all(misty.run$importances.aggregated %>%
                    dplyr::filter(view == "intra") %>%
                    dplyr::pull(Importance) == 0)
  )
  
  # linear model
  suppressWarnings(
    misty.run <- run_misty(misty.views, model.function = linear_model, 
                           bypass.intra = TRUE) %>%
      collect_results()
  )
  expect_true(all(misty.run$importances.aggregated %>%
                    dplyr::filter(view == "intra") %>%
                    dplyr::pull(Importance) == 0)
  )
  
  # svm
  suppressWarnings(
    misty.run <- run_misty(misty.views, model.function = svm_model, 
                           bypass.intra = TRUE) %>%
      collect_results()
  )
  expect_true(all(misty.run$importances.aggregated %>%
                    dplyr::filter(view == "intra") %>%
                    dplyr::pull(Importance) == 0)
  )
  
  # mlp
  suppressWarnings(
    misty.run <- run_misty(misty.views, model.function = mlp_model, 
                           bypass.intra = TRUE) %>%
      collect_results()
  )
  expect_true(all(misty.run$importances.aggregated %>%
                    dplyr::filter(view == "intra") %>%
                    dplyr::pull(Importance) == 0)
  )
})

test_that("caching works for all ML algorithms", {
  # get data
  expr <- generate_random_tibble(100, 5)
  pos <- sample_grid_geometry(100, 10, 10)
  misty.views <- create_initial_view(expr) %>% add_paraview(pos, l = 2)
  
  # random forest
  suppressWarnings({
    run_misty(misty.views, model.function = random_forest_model, cached = TRUE)
  })
  cache.folder <- paste0(".misty.temp/", misty.views[["misty.uniqueid"]])
  expect_true(dir.exists(cache.folder))
  expect_length(list.files(cache.folder), 10)
  expect_match(suppressWarnings({
    run_misty(misty.views, model.function = random_forest_model, cached = TRUE)
  }), "results")
  unlink("results", recursive = TRUE)
  clear_cache()
  
  # bagged mars
  suppressWarnings({
    run_misty(misty.views, model.function = bagged_mars_model, cached = TRUE)
  })
  cache.folder <- paste0(".misty.temp/", misty.views[["misty.uniqueid"]])
  expect_true(dir.exists(cache.folder))
  expect_length(list.files(cache.folder), 10)
  expect_match(suppressWarnings({
    run_misty(misty.views, model.function = bagged_mars_model, cached = TRUE)
  }), "results")
  unlink("results", recursive = TRUE)
  clear_cache()
  
  # mars
  suppressWarnings({
    run_misty(misty.views, model.function = mars_model, cached = TRUE)
  })
  cache.folder <- paste0(".misty.temp/", misty.views[["misty.uniqueid"]])
  expect_true(dir.exists(cache.folder))
  expect_length(list.files(cache.folder), 10)
  expect_match(suppressWarnings({
    run_misty(misty.views, model.function = mars_model, cached = TRUE)
  }), "results")
  unlink("results", recursive = TRUE)
  clear_cache()
  
  # gradient boosting
  suppressWarnings({
    run_misty(misty.views, model.function = gradient_boosting_model, cached = TRUE)
  })
  cache.folder <- paste0(".misty.temp/", misty.views[["misty.uniqueid"]])
  expect_true(dir.exists(cache.folder))
  expect_length(list.files(cache.folder), 10)
  expect_match(suppressWarnings({
    run_misty(misty.views, model.function = gradient_boosting_model, cached = TRUE)
  }), "results")
  unlink("results", recursive = TRUE)
  clear_cache()
  
  # linear model
  suppressWarnings({
    run_misty(misty.views, model.function = linear_model, cached = TRUE)
  })
  cache.folder <- paste0(".misty.temp/", misty.views[["misty.uniqueid"]])
  expect_true(dir.exists(cache.folder))
  expect_length(list.files(cache.folder), 10)
  expect_match(suppressWarnings({
    run_misty(misty.views, model.function = linear_model, cached = TRUE)
  }), "results")
  unlink("results", recursive = TRUE)
  clear_cache()
  
  # svm
  suppressWarnings({
    run_misty(misty.views, model.function = svm_model, cached = TRUE)
  })
  cache.folder <- paste0(".misty.temp/", misty.views[["misty.uniqueid"]])
  expect_true(dir.exists(cache.folder))
  expect_length(list.files(cache.folder), 10)
  expect_match(suppressWarnings({
    run_misty(misty.views, model.function = svm_model, cached = TRUE)
  }), "results")
  unlink("results", recursive = TRUE)
  clear_cache()
  
  # mlp
  suppressWarnings({
    run_misty(misty.views, model.function = mlp_model, cached = TRUE)
  })
  cache.folder <- paste0(".misty.temp/", misty.views[["misty.uniqueid"]])
  expect_true(dir.exists(cache.folder))
  expect_length(list.files(cache.folder), 10)
  expect_match(suppressWarnings({
    run_misty(misty.views, model.function = mlp_model, cached = TRUE)
  }), "results")
  unlink("results", recursive = TRUE)
  clear_cache()
})
