# mistyR utility functions
# Copyleft (ɔ) 2020-2021 Jovan Tanevski <jovan.tanevski@uni-heidelberg.de>

#' Aggregate collected results
#'
#' Helper function
#'
#' @param improvements collected improvements.
#' @param contributions collected contributions.
#' @param importances collected importances.
#'
#' @return a list of improvement stats, contribution stats and
#' aggregated importances.
#'
#' @seealso \code{\link{collect_results}()} to collect results.
#'
#' @noRd
aggregate_results <- function(improvements, contributions, importances) {
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
       dplyr::filter(!stringr::str_starts(.data$view, "p\\.") &
                       .data$view != "intercept") %>%
       dplyr::group_by(.data$target, .data$view) %>%
       dplyr::summarise(mean = mean(.data$value), .groups = "drop_last") %>%
       dplyr::mutate(fraction = abs(.data$mean) / sum(abs(.data$mean))) %>%
       dplyr::ungroup()),
    # p values
    (contributions %>%
       dplyr::filter(stringr::str_starts(.data$view, "p\\.") &
                       !stringr::str_detect(.data$view, "intercept")) %>%
       dplyr::group_by(.data$target, .data$view) %>%
       dplyr::mutate(view = stringr::str_remove(.data$view, "^p\\.")) %>%
       dplyr::summarise(
         p.mean = mean(.data$value),
         p.sd = stats::sd(.data$value),
         .groups = "drop"
       )),
    by = c("target", "view")
  )
  
  importances.aggregated <- importances %>%
    tidyr::unite("PT", "Predictor", "Target") %>%
    dplyr::group_by(.data$view, .data$PT) %>%
    dplyr::summarise(
      Importance = mean(.data$Importance),
      nsamples = dplyr::n(), .groups = "drop"
    ) %>%
    tidyr::separate("PT", c("Predictor", "Target"), sep = "[^(\\.|[:alnum:])]+")
  
  return(list(
    improvements.stats = improvements.stats,
    contributions.stats = contributions.stats,
    importances.aggregated = importances.aggregated
  ))
}


#' Collect and aggregate results
#'
#' Collect and aggregate performance, contribution and importance estimations
#' of a set of raw results produced by \code{\link{run_misty}()}.
#'
#' @param folders Paths to folders containing the raw results from
#'     \code{\link{run_misty}()}.
#'
#' @return List of collected performance, contributions and importances per sample,
#'     performance and contribution statistics and aggregated importances.
#'     \describe{
#'         \item{\var{improvements}}{Long format \code{tibble} with measurements
#'             of performance for each \var{target} and each \var{sample}.
#'             Available performance measures are RMSE and variance explained
#'             (R2) for a model containing only an intrinsic view
#'             (\var{intra.RMSE}, \var{intra.R2}), model with all views
#'             (\var{multi.RMSE}, \var{multi.R2}), gain of RMSE and gain of
#'             variance explained of multi-view model over the intrisic model
#'             where \var{gain.RMSE} is the relative decrease of RMSE in percent,
#'             while \var{gain.R2} is the absolute increase of variance explained
#'             in percent. Each \var{value} represents the mean performance across
#'             folds (k-fold cross-validation). The p values of a one sided
#'             t-test of improvement of performance (\var{p.RMSE}, \var{p.R2})
#'             are also available as a measure.}
#'         \item{\var{improvements.stats}}{Long format \code{tibble} with summary
#'             statistics (mean, standard deviation and coefficient of variation)
#'             for all performance measures for each {target} over all samples.}
#'         \item{\var{contributions}}{Long format \code{tibble} with the values
#'             of the coefficients for each \var{view} in the meta-model, for each
#'             \var{target} and each \var{sample}. The p values for the coefficient
#'             for each view, under the null hypothesis of zero contribution to the
#'             meta model are also available.}
#'         \item{\var{contributions.stats}}{Long format \code{tibble} with summary
#'             statistics for all views per target over all samples. Including
#'             mean coffecient value, fraction of contribution, mean and standard
#'             deviation of p values.}
#'         \item{\var{importances}}{List of view-specific predictor-target
#'         importance tables per sample. The importances in each table are
#'         standardized per target and weighted by the quantile of the coefficient
#'         for the target in that view. Columns other than \var{Predictor}
#'         represent target markers.}
#'         \item{\var{importances.aggregated}}{A list of aggregated view-specific
#'         predictor-target importance tables . Aggregation is
#'         reducing by mean over all samples.}
#'     }
#'
#' @seealso \code{\link{run_misty}()} to train models and
#'     generate results.
#'
#' @examples
#' # Train and collect results for 3 samples in synthetic
#'
#' library(dplyr)
#' library(purrr)
#'
#' data("synthetic")
#'
#' misty.results <- synthetic[seq_len(3)] %>%
#'   imap_chr(~ create_initial_view(.x %>% select(-c(row, col, type))) %>%
#'     add_paraview(.x %>% select(row, col), l = 10) %>%
#'     run_misty(paste0("results/", .y))) %>%
#'   collect_results()
#' str(misty.results)
#' @export
collect_results <- function(folders) {
  samples <- R.utils::getAbsolutePath(folders)
  
  message("\nCollecting improvements")
  improvements <- samples %>%
    furrr::future_map_dfr(function(sample) {
      performance <- readr::read_table(paste0(sample, .Platform$file.sep, "performance.txt"),
                                       na = c("", "NA", "NaN"), col_types = readr::cols()
      ) %>% dplyr::distinct()
      
      performance %>%
        dplyr::mutate(
          sample = sample,
          gain.RMSE = 100 * (.data$intra.RMSE - .data$multi.RMSE) / .data$intra.RMSE,
          gain.R2 = .data$multi.R2 - .data$intra.R2
        )
    }, .progress = TRUE) %>%
    tidyr::pivot_longer(-c(.data$sample, .data$target), names_to = "measure")
  
  
  message("\nCollecting contributions")
  contributions <- samples %>% furrr::future_map_dfr(function(sample) {
    coefficients <- readr::read_table(paste0(sample, .Platform$file.sep, "coefficients.txt"),
                                      na = c("", "NA", "NaN"), col_types = readr::cols()
    ) %>% dplyr::distinct()
    
    coefficients %>%
      dplyr::mutate(sample = sample, .after = "target") %>%
      tidyr::pivot_longer(cols = -c(.data$sample, .data$target), names_to = "view")
  }, .progress = TRUE)
  
  message("\nCollecting importances")
  importances <- samples %>%
    furrr::future_map_dfr(function(sample) {
      targets <- contributions %>%
        dplyr::filter(.data$sample == !!sample) %>%
        dplyr::pull(.data$target) %>%
        unique() %>%
        sort()
      views <- contributions %>%
        dplyr::pull(.data$view) %>%
        unique() %>%
        stringr::str_subset("^p\\.", negate = TRUE) %>%
        stringr::str_subset("^intercept$", negate = TRUE)
      
      # one heatmap per view
      maps <- views %>%
        furrr::future_map_dfr(function(view) {
          all.importances <- targets %>% purrr::map(~ readr::read_csv(paste0(
            sample, .Platform$file.sep, "importances_", .x, "_", view, ".txt"
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
            dplyr::filter(.data$sample == !!sample, view == paste0("p.", !!view)) %>%
            dplyr::mutate(value = 1 - .data$value)
          
          # importances are standardized for each target
          # and multiplied by 1-pval(view)
          all.importances %>%
            purrr::imap_dfc(~
              tibble::tibble(feature = features, zero.imp = 0) %>%
              dplyr::left_join(.x, by = "feature") %>%
              dplyr::arrange(.data$feature) %>%
              dplyr::mutate(
                imp = scale(.data$imp)[, 1],
                !!targets[.y] := .data$zero.imp + (.data$imp *
                   (pvalues %>%
                      dplyr::filter(target == targets[.y]) %>%
                      dplyr::pull(.data$value)))
              )
            %>%
              dplyr::select(targets[.y])) %>%
            dplyr::mutate(Predictor = features) %>%
            tidyr::pivot_longer(
              names_to = "Target",
              values_to = "Importance",
              -.data$Predictor
            ) %>%
            dplyr::mutate(Importance = replace(
              .data$Importance,
              is.nan(.data$Importance), 0
            )) %>%
            dplyr::mutate(view = view, .before = 1)
        }) %>%
        dplyr::mutate(sample = sample, .before = 1)
    }, .progress = TRUE)
  
  message("\nAggregating")
  
  misty.results <- c(
    list(
      improvements = improvements,
      contributions = contributions,
      importances = importances
    ),
    aggregate_results(improvements, contributions, importances)
  )
  
  return(misty.results)
}

#' Aggregate a subset of results
#'
#' @inheritParams collect_results
#'
#' @param misty.results a results list generated by
#'     \code{\link{collect_results}()}.
#'
#' @return the \code{misty.results} list with an added list item
#'     \code{importances.aggregated.subset} containing the aggregated importances
#'     for the subset of \code{folders}.
#'
#' @seealso \code{\link{collect_results}()} to generate a
#'     results list from raw results.
#'
#' @noRd
aggregate_results_subset <- function(misty.results, folders) {
  assertthat::assert_that(("importances" %in% names(misty.results)),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )
  
  normalized.folders <- R.utils::getAbsolutePath(folders)
  # check if folders are in names of misty.results
  assertthat::assert_that(all(normalized.folders %in%
    (misty.results$importances %>% dplyr::pull(.data$sample))),
    msg = "The provided results list doesn't contain information about some of
    the requested result folders. Consider using collect_results()."
  )
  
  message("Aggregating subset")
  importances.aggregated.subset <- misty.results$importances %>%
    dplyr::filter(.data$sample %in% normalized.folders) %>%
    tidyr::unite("PT", "Predictor", "Target") %>%
    dplyr::group_by(.data$view, .data$PT) %>%
    dplyr::summarise(
      Importance = mean(.data$Importance),
      nsamples = dplyr::n(), .groups = "drop"
    ) %>%
    tidyr::separate("PT", c("Predictor", "Target"))
  
  misty.results[["importances.aggregated.subset"]] <- importances.aggregated.subset
  
  return(misty.results)
}

#' Clear cached objects
#'
#' Purge the cache or clear the cached objects for a single sample.
#'
#' The cached objects are removed from disk and cannot be retrieved. Whenever
#' possible specifying an \code{id} is reccomended. If \code{id = NULL} all
#' contents of the folder \file{.misty.temp} will be removed.
#'
#' @param id the unique id of the sample.
#'
#' @return None (\code{NULL})
#'
#' @examples
#' clear_cache("b98ad35f4e671871cba35f2155228612")
#'
#' clear_cache()
#' @export
clear_cache <- function(id = NULL) {
  cache.folder <- R.utils::getAbsolutePath(".misty.temp")
  if (is.null(id)) {
    if (dir.exists(cache.folder)) {
      unlink(cache.folder, recursive = TRUE)
    } else {
      warning("Cache folder doesn't exist.")
    }
  } else {
    sample.cache.folder <- paste0(cache.folder, .Platform$file.sep, id)
    if (dir.exists(sample.cache.folder)) {
      unlink(sample.cache.folder, recursive = TRUE)
    } else {
      warning("Cache folder for requested id doesn't exist.")
    }
  }
}

#' Removes empty cache folders.
#'
#' @return None (\code{NULL})
#'
#' @noRd
sweep_cache <- function() {
  cache.folder <- R.utils::getAbsolutePath(".misty.temp")
  if (dir.exists(cache.folder)) {
    list.files(cache.folder, full.names = TRUE) %>%
      purrr::walk(function(path) {
        if (length(list.files(path)) == 0) {
          unlink(path, recursive = TRUE)
        }
      })
    
    if (length(list.files(cache.folder, full.names = TRUE)) == 0) {
      clear_cache()
    }
  }
}


#' Extract signatures from the results
#'
#' Signature is a representation of each sample in the space of mistyR results.
#'
#' The performance signature of each sample is a concatenation of the estimated
#' values of variance explained using only the intraview, the variance explained
#' by the multiview model and the gain in variance explained for each marker.
#' The performance signature vector for each sample available in
#' \code{misty.results} is of length \eqn{\textrm{markers} \cdot 3}{markers x 3}.
#'
#' The contribution signature of each sample is a concatenation of the estimated
#' fraction of contribution of each view for each marker.
#' The contribution signature vector for each sample available in
#' \code{misty.results} is of length
#' \eqn{\textrm{markers} \cdot \textrm{views}}{markers x views}.
#'
#' The importance signature of each sample is a concatenation of the estimated
#' and weighted importances for each predictor-target marker pair from all views.
#' The importance signature vector for each sample available in
#' \code{misty.results} is of length
#' \eqn{\textrm{markers}^2 \cdot \textrm{views}}{markers^2 x views}.
#'
#' @inheritParams plot_interaction_heatmap
#' @param type type of signature to extract from the results.
#'
#' @return A table with one row per sample from \code{misty.results} representing
#' its signature.
#'
#' @seealso \code{\link{collect_results}()} to generate a
#'     results list from raw results.
#'
#' @examples
#' library(dplyr)
#'
#' misty.results <-
#'   list.files("results", full.names = TRUE) %>% collect_results()
#'
#' extract_signature(misty.results, "performance")
#' @export
extract_signature <- function(misty.results,
                              type = c(
                                "performance", "contribution", "importance"
                              ), trim = -Inf, trim.measure = c(
                                "gain.R2", "multi.R2", "intra.R2",
                                "gain.RMSE", "multi.RMSE", "intra.RMSE"
                              )) {
  signature.type <- match.arg(type)
  
  trim.measure.type <- match.arg(trim.measure)
  
  assertthat::assert_that(
    all(c(
      "improvements", "contributions",
      "importances", "importances.aggregated"
    ) %in%
      names(misty.results)),
    msg = "The provided result list is malformed.
    Consider using collect_results()."
  )
  
  inv <- sign((stringr::str_detect(trim.measure.type, "gain") |
    stringr::str_detect(trim.measure.type, "RMSE", negate = TRUE)) - 0.5)
  
  targets <- misty.results$improvements.stats %>%
    dplyr::filter(
      .data$measure == trim.measure.type,
      inv * .data$mean >= inv * trim
    ) %>%
    dplyr::pull(.data$target)
  
  switch(signature.type,
         "performance" = {
           target.intersection <- misty.results$improvements %>%
             dplyr::group_by(.data$sample) %>%
             dplyr::summarize(ts = list(unique(.data$target))) %>%
             dplyr::pull(.data$ts) %>%
             purrr::reduce(intersect) %>%
             intersect(targets)
           
           misty.results$improvements %>%
             dplyr::filter(
               .data$target %in% target.intersection,
               stringr::str_ends(.data$measure, "R2"),
               !stringr::str_ends(.data$measure, "p.R2")
             ) %>%
             tidyr::unite("Feature", .data$target, .data$measure) %>%
             # grouping necessary?
             dplyr::group_by(.data$sample) %>%
             tidyr::pivot_wider(names_from = "Feature", values_from = "value") %>%
             dplyr::ungroup()
         },
         "contribution" = {
           target.intersection <- misty.results$contributions %>%
             dplyr::group_by(.data$sample) %>%
             dplyr::summarize(ts = list(unique(.data$target))) %>%
             dplyr::pull(.data$ts) %>%
             purrr::reduce(intersect) %>%
             intersect(targets)
           
           misty.results$contributions %>%
             dplyr::filter(
               .data$target %in% target.intersection,
               !stringr::str_starts(.data$view, "p\\."),
               !stringr::str_detect(.data$view, "intercept")
             ) %>%
             dplyr::group_by(.data$sample, .data$target) %>%
             dplyr::mutate(frac = abs(.data$value) / sum(abs(.data$value)), value = NULL) %>%
             tidyr::unite("Feature", .data$view, .data$target) %>%
             tidyr::pivot_wider(names_from = "Feature", values_from = "frac") %>%
             dplyr::ungroup()
         },
         "importance" = {
           views <- misty.results$importances %>%
             dplyr::pull(.data$view) %>%
             unique()
           
           views %>%
             purrr::map(function(view) {
               view.importances <- misty.results$importances %>%
                 dplyr::filter(.data$view == !!view, !is.na(.data$Importance))
               target.intersection <- view.importances %>%
                 dplyr::group_by(.data$sample) %>%
                 dplyr::summarize(ts = list(unique(.data$Target))) %>%
                 dplyr::pull(.data$ts) %>%
                 purrr::reduce(intersect) %>%
                 intersect(targets)
               predictor.intersection <- view.importances %>%
                 dplyr::group_by(.data$sample) %>%
                 dplyr::summarize(ts = list(unique(.data$Predictor))) %>%
                 dplyr::pull(.data$ts) %>%
                 purrr::reduce(intersect)
               
               view.importances %>%
                 dplyr::filter(
                   .data$Predictor %in% predictor.intersection,
                   .data$Target %in% target.intersection
                 ) %>%
                 tidyr::unite("vPT", .data$view, .data$Predictor, .data$Target) %>%
                 tidyr::pivot_wider(names_from = "vPT", values_from = "Importance")
             }) %>%
             purrr::reduce(dplyr::full_join, by = "sample")
         }
  )
}
