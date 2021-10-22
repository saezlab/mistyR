test_that("rename_views() behaves correctly", {
  expr <- generate_random_tibble(30, 5)
  pos <- sample_grid_geometry(30, 10, 10)
  misty.views <- create_initial_view(expr) %>%
    add_paraview(positions = pos, l = 10, family = "gaussian")
  tname <- "test.view"
  expect_false(tname %in% names(misty.views))
  misty.views <- misty.views %>% rename_view("paraview.10", tname)
  expect_true(tname %in% names(misty.views))
  expect_false("paraview.10" %in% names(misty.views))
  expect_error(misty.views %>% rename_view("bar", "foo"))
})

test_that("filter_views() behaves correctly", {
  expr <- generate_random_tibble(30, 5) %>%
    dplyr::mutate(expr6 = c(rep(5, 15), rep(10, 15)))
  pos <- sample_grid_geometry(30, 10, 10)
  misty.views <- create_initial_view(expr) %>%
    add_paraview(positions = pos, l = 10, family = "gaussian")
  truncated_views <- misty.views %>%
    filter_views(rows = 1:10)
  expect_equal(nrow(truncated_views$intraview$data), 10)
  expect_equal(nrow(truncated_views$paraview.10$data), 10)
  truncated_views <- misty.views %>%
    filter_views(NA, "intraview", expr6 < 8)
  expect_equal(nrow(truncated_views$intraview$data), 15)
  expect_equal(nrow(truncated_views$paraview.10$data), 15)
  expect_error(misty.views %>% filter_views(NA, view = "bar"))
  expect_error(misty.views %>% filter_views(NA, view = "intraview", expr7 < 8))
})

test_that("select_markers() behaves correctly", {
  expr <- tibble::tibble(
    prod1 = rnorm(30, 10, 2),
    prod2 = rnorm(30, 20, 2),
    lig1 = runif(30, 1, 5),
    lig2 = runif(30, 3, 7)
  )
  misty.views <- create_initial_view(expr)
  truncated_views <- misty.views %>%
    select_markers("intraview", matches("^prod[0-9]+$"))
  expect_equal(ncol(truncated_views$intraview$data), 2)
  truncated_views <- misty.views %>%
    select_markers("intraview", matches("^[a-z]+1$"))
  expect_equal(ncol(truncated_views$intraview$data), 2)
  expect_error(misty.views %>% select_markers("bar", startswith("prod")))
})
