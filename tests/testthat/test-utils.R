test_that("clean_cache cleans correctly", {
  expr <- generate_random_tibble(30, 5)
  pos <- sample_grid_geometry(30, 10, 10)
  misty.views <- create_initial_view(expr) %>%
    add_paraview(pos, 2, cached = TRUE)
  expect_true(dir.exists(".misty.temp"))
  expect_length(list.files(".misty.temp", ), 1)
  clear_cache(misty.views[["misty.uniqueid"]])
  expect_length(list.files(".misty.temp", ), 0)
  clear_cache()
  expect_false(dir.exists(".misty.temp"))
  expect_warning(clear_cache("nonexistent"))
  expect_warning(clear_cache())
})

test_that("sweep_cache sweeps empty folders only", {
  dir.create(".misty.temp")
  sweep_cache()
  expect_false(dir.exists(".misty.temp"))
  dir.create(".misty.temp/emptyfolder", recursive = TRUE)
  dir.create(".misty.temp/nonemptyfolder", recursive = TRUE)
  file.create(".misty.temp/nonemptyfolder/file")
  sweep_cache()
  expect_length(list.files(".misty.temp"), 1)
  clear_cache()
})

test_that("collect_results creates expected structure", {
  targets <- 5
  samples <- 3
  suppressWarnings({
    seq_len(samples) %>% walk(function(id) {
      expr <- generate_random_tibble(30, targets, id)
      pos <- sample_grid_geometry(30, 10, 10, id)
      create_initial_view(expr) %>%
        add_paraview(pos, 2) %>%
        run_misty(paste0("results/results", id))
    })
  })
  misty.results <- collect_results(list.files("results", full.names = TRUE))
  expect_length(misty.results, 6)
  expect_equal(nrow(misty.results$improvements), samples*targets*8)
  expect_equal(nrow(misty.results$improvements.stats), targets*6)
  expect_equal(nrow(misty.results$contributions), samples*targets*6)
  expect_equal(nrow(misty.results$contributions.stats), targets*2)
  expect_equal(nrow(misty.results$importances), samples*2*targets^2)
  expect_equal(nrow(misty.results$importances.aggregated), 2*targets^2)
  subset.results <- aggregate_results_subset(misty.results, 
                                             c("results/results1",
                                               "results/results2"))
  expect_length(subset.results, 7)
  expect_equal(nrow(subset.results$importances.aggregated.subset), 2*targets^2)
  unlink("results", recursive = TRUE)
})
