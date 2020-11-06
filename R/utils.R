#' Clear cached models
#'
#' @param singleid
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
clear_cache <- function(singleid = NULL) {
  if (is.null(singleid)) {
    if (unlink(".misty.temp", recursive = TRUE) == 0) {
      message("Cache cleared")
    } else {
      message("Failed to clear cache")
    }
  } else {
    if (unlink(paste0(
      ".misty.temp", .Platform$file.sep, singleid
    ), recursive = TRUE) == 0) {
      message("Cache cleared\n")
    } else {
      message("Failed to clear cache\n")
    }
  }
}


#' Collect and aggregate MISTy results
#'
#' @param folders
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
collect_results <- function(folders) {
  images <- folders[dir.exists(folders)]

  message(paste("Collecting results from", paste(images, collapse = " ")))

  message("\nCollecting improvements")
  improvements <- images %>%
    furrr::future_map_dfr(function(image) {
      performance <- readr::read_delim(paste0(image, .Platform$file.sep, "performance.txt"),
        delim = " ", col_types = readr::cols()
      ) %>% dplyr::distinct()

      performance %>%
        dplyr::mutate(
          image = image,
          gain.RMSE = 100 * (.data$intra.RMSE - .data$multi.RMSE) / .data$intra.RMSE,
          gain.R2 = 100 * (.data$multi.R2 - .data$intra.R2),
        )
    }, .progress = TRUE) %>%
    tidyr::pivot_longer(-c(.data$image, .data$target), names_to = "measure")


  message("\nCollecting contributions")
  contributions <- images %>% furrr::future_map_dfr(function(image) {
    coefficients <- readr::read_delim(paste0(image, .Platform$file.sep, "coefficients.txt"),
      delim = " ", col_types = readr::cols()
    ) %>% dplyr::distinct()

    coefficients %>%
      dplyr::mutate(image = image, .after = "target") %>%
      tidyr::pivot_longer(cols = -c(.data$image, .data$target), names_to = "view")
  }, .progress = TRUE)

  improvements.stats <- improvements %>%
    dplyr::filter(!stringr::str_starts(.data$measure, "p\\.")) %>%
    dplyr::group_by(.data$target, .data$measure) %>%
    dplyr::summarise(
      mean = mean(.data$value), sd = stats::sd(.data$value),
      cv = .data$sd / .data$mean, .groups = "drop"
    )


  contributions.stats <- dplyr::inner_join(
    # mean coefficients
    (contributions %>%
      dplyr::filter(!stringr::str_starts(.data$view, "p\\.") & .data$view != "intercept") %>%
      dplyr::group_by(.data$target, .data$view) %>%
      dplyr::summarise(mean = mean(.data$value), .groups = "drop_last") %>%
      dplyr::mutate(fraction = abs(.data$mean) / sum(abs(.data$mean))) %>%
      dplyr::ungroup()),
    # p values
    (contributions %>%
      dplyr::filter(stringr::str_starts(.data$view, "p\\.") & !stringr::str_detect(.data$view, "intercept")) %>%
      dplyr::group_by(.data$target, .data$view) %>%
      dplyr::mutate(view = stringr::str_remove(.data$view, "^p\\.")) %>%
      dplyr::summarise(p.mean = mean(.data$value), p.sd = stats::sd(.data$value), .groups = "drop")),
    by = c("target", "view")
  )

  message("\nCollecting importances")
  importances <- images %>%
    furrr::future_map(function(image) {
      targets <- contributions.stats %>%
        dplyr::pull(.data$target) %>%
        unique() %>%
        sort()
      views <- contributions.stats %>%
        dplyr::pull(.data$view) %>%
        unique()

      # one heatmap per view
      maps <- views %>%
        furrr::future_map(function(view) {
          all.importances <- targets %>% purrr::map(~ readr::read_csv(paste0(
            image, .Platform$file.sep, "importances_", .x, "_", view, ".txt"
          ),
          col_types = readr::cols()
          ) %>%
            dplyr::distinct() %>%
            dplyr::rename(feature = target))

          features <- all.importances %>%
            purrr::map(~ .x$feature) %>%
            unlist() %>%
            unique() %>%
            sort()

          pvalues <- contributions %>%
            dplyr::filter(image == !!image, view == paste0("p.", !!view)) %>%
            dplyr::mutate(value = 1 - .data$value)

          # importances are standardized for each target an multiplied by 1-pval(view)
          all.importances %>%
            purrr::imap_dfc(~
            tibble::tibble(feature = features, zero.imp = 0) %>%
              dplyr::left_join(.x, by = "feature") %>%
              dplyr::arrange(.data$feature) %>%
              dplyr::mutate(imp = scale(.data$imp)[, 1], !!targets[.y] := .data$zero.imp + (.data$imp *
                (pvalues %>%
                  dplyr::filter(target == targets[.y]) %>%
                  dplyr::pull(.data$value))))
              %>%
              dplyr::select(targets[.y])) %>%
            dplyr::mutate(Predictor = features)
        }) %>%
        `names<-`(views)
    }, .progress = TRUE) %>%
    `names<-`(images)

  message("\nAggregating")
  importances.aggregated <- importances %>%
    purrr::reduce(function(acc, l) {
      acc %>% purrr::map2(l, ~ (((.x %>% dplyr::select(-.data$Predictor)) +
        (.y %>% dplyr::select(-.data$Predictor))) %>%
        dplyr::mutate(Predictor = .x %>% dplyr::pull(.data$Predictor))))
    }) %>%
    purrr::map(~ .x %>% dplyr::mutate_if(is.numeric, ~ . / length(images)))

  return(list(
    improvements = improvements,
    improvements.stats = improvements.stats,
    contributions = contributions,
    contributions.stats = contributions.stats,
    importances = importances,
    importances.aggregated = importances.aggregated
  ))
}


#' Title
#'
#' @param misty.results
#' @param folders
#'
#' @return
#' @export
#'
#' @examples
#' # TBD
aggregate_results_subset <- function(misty.results, folders) {
  assertthat::assert_that(("importances" %in% names(misty.results)),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )

  # check if folders are in names of misty.results
  assertthat::assert_that(all(folders %in% names(misty.results$importances)),
    msg = "The provided result list doesn't contain information about some of the requested result folders.
    Consider using collect_results()."
  )

  message("Aggregating subset")
  importances.aggregated.subset <- rlist::list.subset(misty.results$importances, folders) %>%
    purrr::reduce(function(acc, l) {
      acc %>% purrr::map2(l, ~ (((.x %>% dplyr::select(-Predictor)) + (.y %>% dplyr::select(-Predictor))) %>%
        dplyr::mutate(Predictor = .x %>% dplyr::pull(Predictor))))
    }) %>%
    purrr::map(~ .x %>% dplyr::mutate_if(is.numeric, ~ . / length(folders)))

  misty.results[["importances.aggregated.subset"]] <- importances.aggregated.subset

  return(misty.results)
}
