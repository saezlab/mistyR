# misty.results is a list obtained by running the function collect_results

#' Plot observed performance and improvement per target
#'
#' @param misty.results
#' @param measure
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
plot_improvement_stats <- function(misty.results, measure = "gain.R2") {
  assertthat::assert_that(("improvements.stats" %in% names(misty.results)),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )

  plot.data <- misty.results$improvements.stats %>%
    dplyr::filter(measure == !!measure)

  assertthat::assert_that(nrow(plot.data) > 0,
    msg = "The selected measure cannot be found in the results table."
  )

  set2.orange <- "#FC8D62"

  results.plot <- ggplot2::ggplot(
    plot.data,
    ggplot2::aes(
      x = stats::reorder(.data$target, -.data$mean),
      y = .data$mean
    )
  ) +
    ggplot2::geom_pointrange(ggplot2::aes(
      ymin = .data$mean - .data$sd,
      ymax = .data$mean + .data$sd
    )) +
    ggplot2::geom_point(color = set2.orange) +
    ggplot2::theme_classic() +
    ggplot2::ylab(measure) +
    ggplot2::xlab("Target") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  print(results.plot)

  invisible(misty.results)
}

#' Title
#'
#' @param misty.results
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
plot_view_contributions <- function(misty.results) {
  assertthat::assert_that(("contributions.stats" %in% names(misty.results)),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )

  plot.data <- misty.results$contributions.stats


  results.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = .data$target, y = .data$fraction)) +
    ggplot2::geom_col(ggplot2::aes(group = .data$view, fill = .data$view)) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::ylab("Contribution") +
    ggplot2::xlab("Target") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  print(results.plot)

  invisible(misty.results)
}

#' Title
#'
#' @param misty.results
#' @param view
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
plot_interaction_heatmap <- function(misty.results, view, cutoff = 1) {
  assertthat::assert_that(("importances.aggregated" %in% names(misty.results)),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )

  assertthat::assert_that((view %in% names(misty.results$importances.aggregated)),
    msg = "The selected view cannot be found in the results table."
  )

  plot.data <- misty.results$importances.aggregated[[view]] %>%
    tidyr::pivot_longer(names_to = "Target", values_to = "Importance", -.data$Predictor)

  set2.blue <- "#8DA0CB"

  results.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = .data$Predictor, y = .data$Target)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$Importance)) +
    ggplot2::scale_fill_gradient2(low = "white", mid = "white", high = set2.blue, midpoint = cutoff) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::ggtitle(view)

  print(results.plot)

  invisible(misty.results)
}

#' Title
#'
#' @param misty.results
#' @param from.view
#' @param to.view
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
plot_contrast_heatmap <- function(misty.results, from.view, to.view, cutoff = 1) {
  assertthat::assert_that(("importances.aggregated" %in% names(misty.results)),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )

  assertthat::assert_that((from.view %in% names(misty.results$importances.aggregated)),
    msg = "The selected from.view cannot be found in the results table."
  )

  assertthat::assert_that((to.view %in% names(misty.results$importances.aggregated)),
    msg = "The selected to.view cannot be found in the results table."
  )

  mask <- ((misty.results$importances.aggregated[[from.view]] %>% dplyr::select(-.data$Predictor)) < cutoff) &
    ((misty.results$importances.aggregated[[to.view]] %>% dplyr::select(-.data$Predictor)) >= cutoff)

  masked <- ((misty.results$importances.aggregated[[to.view]] %>%
    tibble::column_to_rownames("Predictor")) * mask)

  plot.data <- masked %>%
    dplyr::slice(which(masked %>% rowSums(na.rm = TRUE) > 0)) %>%
    dplyr::select(which(masked %>% colSums(na.rm = TRUE) > 0)) %>%
    tibble::rownames_to_column("Predictor") %>%
    tidyr::pivot_longer(names_to = "Target", values_to = "Importance", -.data$Predictor)

  set2.blue <- "#8DA0CB"

  results.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = .data$Predictor, y = .data$Target)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$Importance)) +
    ggplot2::scale_fill_gradient2(low = "white", mid = "white", high = set2.blue, midpoint = cutoff) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::ggtitle(paste0(to.view, " - ", from.view))

  print(results.plot)

  invisible(misty.results)
}

#' Title
#'
#' @param misty.results
#' @param view
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
plot_interaction_communities <- function(misty.results, view, cutoff = 1) {
  assertthat::assert_that(("importances.aggregated" %in% names(misty.results)),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )

  assertthat::assert_that((view %in% names(misty.results$importances.aggregated)),
    msg = "The selected view cannot be found in the results table."
  )

  assertthat::assert_that(
    all(misty.results$importances.aggregated[[view]] %>% dplyr::select(-.data$Predictor) %>% colnames() ==
      misty.results$importances.aggregated[[view]] %>% dplyr::pull(.data$Predictor)),
    msg = "The predictor and target markers in the view must match."
  )

  assertthat::assert_that(requireNamespace("igraph", quietly = TRUE),
    msg = "The package igraph is required to calculate the interaction communities."
  )

  A <- misty.results$importances.aggregated[[view]] %>%
    dplyr::select(-.data$Predictor) %>%
    as.matrix()
  A[A < cutoff | is.na(A)] <- 0

  # Workaround: fix binding check for variable . in this function context.
  # The reference to . below works just fine in its own scope
  . <- NULL

  G <- igraph::graph.adjacency(A, mode = "plus", weighted = TRUE) %>%
    igraph::set.vertex.attribute("name", value = names(igraph::V(.)))

  C <- igraph::cluster_louvain(G)

  layout <- igraph::layout.fruchterman.reingold(G)
  igraph::plot.igraph(G,
    layout = layout, mark.groups = C, main = view, vertex.size = 4,
    vertex.color = "black", vertex.label.dist = 1
  )

  invisible(misty.results)
}


#' Title
#'
#' @param misty.results.from
#' @param misty.results.to
#' @param view
#' @param cutoff.from
#' @param cutoff.to
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
plot_contrast_results <- function(misty.results.from, misty.results.to, views = NULL, cutoff.from = 1, cutoff.to = 1) {
  assertthat::assert_that(("importances.aggregated" %in% names(misty.results.from)),
    msg = "The first provided result list is malformed. Consider using collect_results()."
  )

  assertthat::assert_that(("importances.aggregated" %in% names(misty.results.to)),
    msg = "The second provided result list is malformed. Consider using collect_results()."
  )

  if (is.null(views)) {
    assertthat::assert_that(rlang::is_empty(setdiff(
      names(misty.results.from$importances.aggregated),
      names(misty.results.to$importances.aggregated)
    )),
    msg = "The requested views do not exist in both result lists."
    )
    views <- names(misty.results.from$importances.aggregated)
  } else {
    assertthat::assert_that(all(views %in% names(misty.results.from$importances.aggregated)) &
      all(views %in% names(misty.results.to$importances.aggregated)),
    msg = "The requested views do not exist in both result lists."
    )
  }

  assertthat::assert_that(
    all(views %>% purrr::map_lgl(function(current.view) {
      rlang::is_empty(setdiff(
        misty.results.from$importances.aggregated[[current.view]] %>% colnames(),
        misty.results.to$importances.aggregated[[current.view]] %>% colnames()
      )) &
        rlang::is_empty(setdiff(
          misty.results.from$importances.aggregated[[current.view]] %>% dplyr::pull(.data$Predictor),
          misty.results.to$importances.aggregated[[current.view]] %>% dplyr::pull(.data$Predictor)
        ))
    })),
    msg = "Incompatible predictors and targets."
  )

  views %>% purrr::walk(function(current.view) {
    mask <- ((misty.results.from$importances.aggregated[[current.view]] %>% dplyr::select(-.data$Predictor)) < cutoff.from) &
      ((misty.results.to$importances.aggregated[[current.view]] %>% dplyr::select(-.data$Predictor)) >= cutoff.to)

    masked <- ((misty.results.to$importances.aggregated[[current.view]] %>%
      tibble::column_to_rownames("Predictor")) * mask)

    plot.data <- masked %>%
      dplyr::slice(which(masked %>% rowSums(na.rm = TRUE) > 0)) %>%
      dplyr::select(which(masked %>% colSums(na.rm = TRUE) > 0)) %>%
      tibble::rownames_to_column("Predictor") %>%
      tidyr::pivot_longer(names_to = "Target", values_to = "Importance", -.data$Predictor)

    set2.blue <- "#8DA0CB"

    results.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = .data$Predictor, y = .data$Target)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$Importance)) +
      ggplot2::scale_fill_gradient2(low = "white", mid = "white", high = set2.blue, midpoint = cutoff.to) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::ggtitle(current.view)

    print(results.plot)
  })

  invisible(misty.results.from)
}
