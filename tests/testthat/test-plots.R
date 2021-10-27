expr <- generate_random_tibble(30, 5)
pos <- sample_grid_geometry(30, 10, 10)
suppressWarnings({
  create_initial_view(expr) %>%
    add_paraview(pos, 2) %>%
    run_misty()
})

misty.results <- collect_results("results")

test_that("improvement_stats runs successfully", {
  expect_invisible(suppressWarnings(plot_improvement_stats(misty.results)))
  expect_invisible(suppressWarnings(
    plot_improvement_stats(misty.results, measure = "intra.R2")
  ))
})

test_that("view_contributions runs successfully", {
  expect_invisible(suppressWarnings(plot_view_contributions(misty.results)))
})

test_that("interaction_heatmap runs successfully", {
  expect_invisible(suppressWarnings(plot_interaction_heatmap(
    misty.results,
    "para.2"
  )))
  expect_invisible(suppressWarnings(plot_interaction_heatmap(
    misty.results,
    "intra"
  )))
  expect_invisible(suppressWarnings(plot_interaction_heatmap(
    misty.results, clean=TRUE,
    "intra"
  )))
})

test_that("contrast_heatmap runs successfully", {
  expect_invisible(suppressWarnings(plot_contrast_heatmap(
    misty.results,
    "intra",
    "para.2",
    0.5
  )))
})

test_that("interaction_communities runs successfully", {
  expect_invisible(suppressWarnings(plot_interaction_communities(
    misty.results,
    "para.2",
    0.5
  )))
})

test_that("contrast_results runs successfully", {
  expr <- generate_random_tibble(30, 5, 43)
  pos <- sample_grid_geometry(30, 10, 10, 43)
  suppressWarnings({
    create_initial_view(expr) %>%
      add_paraview(pos, 2) %>%
      run_misty("results2")
  })
  misty.results2 <- collect_results("results2")
  expect_invisible(suppressWarnings(plot_contrast_results(misty.results,
    misty.results2,
    c("intra", "para.2"),
    cutoff.from = 0.5,
    cutoff.to = 0.5
  )))
  expect_invisible(plot_contrast_results(misty.results,
    misty.results2,
    cutoff.from = 0.5,
    cutoff.to = 0.5
  ))
  misty.results2$importances.aggregated <- 
    misty.results2$importances.aggregated %>% dplyr::filter(view != "para.2") 
  expect_error(plot_contrast_results(misty.results,
    misty.results2,
    cutoff.from = 0.5,
    cutoff.to = 0.5
  ))
})

unlink(c("results", "results2"), recursive = TRUE)
if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
