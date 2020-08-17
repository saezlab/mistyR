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
      cat("Cache cleared\n")
    } else {
      cat("Failed to clear cache\n")
    }
  } else {
    if (unlink(paste0(
      ".misty.temp", .Platform$file.sep, singleid
    ), recursive = TRUE) == 0) {
      cat("Cache cleared\n")
    } else {
      cat("Failed to clear cache\n")
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
  
  improvements <- images %>%
    purrr::map_dfr(function(image) {
      performance <- readr::read_delim(paste0(image, .Platform$file.sep, "performance.txt"),
        delim = " ", col_types = readr::cols()
      ) %>% dplyr::distinct()

      performance %>%
        dplyr::mutate(
          image = image,
          gain.RMSE = 100 * (intra.RMSE - multi.RMSE) / intra.RMSE,
          gain.R2 = 100 * (multi.R2 - intra.R2),
        )
    }) %>%
    tidyr::pivot_longer(-c(image, target), names_to = "measure")

  contributions <- images %>% purrr::map_dfr(function(image) {
    coefficients <- readr::read_delim(paste0(image, .Platform$file.sep, "coefficients.txt"),
      delim = " ", col_types = readr::cols()
    ) %>% dplyr::distinct()

    coefficients %>%
      dplyr::mutate(image = image, .after = "target") %>%
      tidyr::pivot_longer(cols = -c(image, target), names_to = "view")
  })

  improvements.stats <- improvements %>%
    dplyr::filter(!stringr::str_starts(measure, "p\\.")) %>%
    dplyr::group_by(target, measure) %>%
    dplyr::summarise(mean = mean(value), sd = sd(value), cv = sd / mean, .groups = "drop")


  contributions.stats <- dplyr::inner_join(
    # mean coefficients
    (contributions %>%
      dplyr::filter(!stringr::str_starts(view, "p\\.") & view != "intercept") %>%
      dplyr::group_by(target, view) %>%
      dplyr::summarise(mean = mean(value), .groups = "drop_last") %>%
      dplyr::mutate(fraction = abs(mean) / sum(abs(mean))) %>%
      dplyr::ungroup()),
    # p values
    (contributions %>%
      dplyr::filter(stringr::str_starts(view, "p\\.") & !stringr::str_detect(view, "intercept")) %>%
      dplyr::group_by(target, view) %>%
      dplyr::mutate(view = stringr::str_remove(view, "^p\\.")) %>%
      dplyr::summarise(p.mean = mean(value), p.sd = sd(value), .groups = "drop")),
    by = c("target", "view")
  )

  importances <- images %>%
    purrr::map(function(image) {
      targets <- contributions.stats %>%
        dplyr::pull(target) %>%
        unique() %>%
        sort()
      views <- contributions.stats %>%
        dplyr::pull(view) %>%
        unique()

      # one heatmap per view
      maps <- views %>%
        purrr::map(function(view) {
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
            dplyr::mutate(value = 1 - value)

          # importances are standardized for each target an multiplied by 1-pval(view)
          all.importances %>%
            purrr::imap_dfc(~
            tibble::tibble(feature = features, zero.imp = 0) %>%
              dplyr::left_join(.x, by = "feature") %>%
              dplyr::arrange(feature) %>%
              dplyr::mutate(imp = scale(imp)[, 1], !!targets[.y] := zero.imp + (imp *
                (pvalues %>%
                  dplyr::filter(target == targets[.y]) %>%
                  dplyr::pull(value))))
              %>%
              dplyr::select(targets[.y])) %>%
            dplyr::mutate(Predictor = features)
        }) %>%
        `names<-`(views)
    }) %>%
    `names<-`(images)

  importances.aggregated <- importances %>%
    purrr::reduce(function(acc, l) {
      acc %>% purrr::map2(l, ~ (((.x %>% dplyr::select(-Predictor)) + (.y %>% dplyr::select(-Predictor))) %>%
        dplyr::mutate(Predictor = .x %>% dplyr::pull(Predictor))))
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
