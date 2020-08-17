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

  results.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = reorder(target, -mean), y = mean)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = mean - sd, ymax = mean + sd)) +
    ggplot2::theme_classic() +
    ggplot2::ylab(measure) +
    ggplot2::xlab("Target") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  print(results.plot)
  
  invisible(misty.results)
}
