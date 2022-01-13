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

test_that("run_misty is reproducible", {
  suppressWarnings({
    hash.results1 <- list.files(run_misty(misty.views, "results1", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results2 <- list.files(run_misty(misty.views, "results2", seed = 1),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
    hash.results3 <- list.files(run_misty(misty.views, "results3", seed = 2),
      full.names = TRUE
    ) %>%
      map_chr(~ digest::digest(.x, file = TRUE))
  })
  expect_equal(hash.results1, hash.results2)
  expect_true(all(hash.results1 != hash.results3))
  unlink(paste0("results", seq_len(3)), recursive = TRUE)
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

test_that("run_misty models are cached and retrieved", {
  suppressWarnings({
    run_misty(misty.views, cached = TRUE)
  })
  cache.folder <- paste0(".misty.temp/", misty.views[["misty.uniqueid"]])
  expect_true(dir.exists(cache.folder))
  expect_length(list.files(cache.folder), 10)
  expect_match(suppressWarnings({
    run_misty(misty.views, cached = TRUE)
  }), "results")
  unlink("results", recursive = TRUE)
  clear_cache()
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
