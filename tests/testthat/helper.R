# generate tibble with normally distributed variables with random means on scale
# [1, 10] and sd equal to the log of the mean
generate_random_tibble <- function(rows, columns, seed = 42) {
  withr::with_seed(seed, {
    means <- runif(columns, 1, 10)
    seq_len(columns) %>%
      purrr::map_dfc(~ tibble::as_tibble_col(
        rnorm(rows, means[.x], log(means[.x])),
        column_name = paste0("expr", .x)
      ))
  })
}

# sample without replacement from a grid.rows * grid.columns grid
sample_grid_geometry <- function(rows, grid.rows, grid.columns, seed = 42) {
  withr::with_seed(seed, {
    row.ind <- sample.int(grid.rows * grid.columns, rows)
    tibble::tibble(
      row = ceiling(row.ind / grid.rows),
      col = row.ind %% grid.rows
    ) %>%
      dplyr::mutate(col = ifelse(col == 0, grid.rows, col))
  })
}
