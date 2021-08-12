# mistyR view manipulation functions
# Copyright (c) 2021 Jovan Tanevski <jovan.tanevski@uni-heidelberg.de>

#' Filter spatial units
#'
#' Select, remove (or duplicate) rows from all views in a composition by their
#' row locations or according to conditions based on a specific view.
#'
#' The values in \code{rows} have priority over the other parameters. If
#' \code{rows} doesn't contain integer values then filtering
#' is performed based on the view specified in \code{view} and expressions
#' (\code{\link[dplyr:dplyr_data_masking]{...}}) returning logical values
#' defined in terms of the variables in \code{view}.
#'
#' @param current.views the current view composition.
#' @param rows row (integer) location; positive values to keep (duplicate)
#' and/or negative to remove.
#' @param view the name of the view to be used for filtering.
#' @param ... logical expressions defined in terms of the variables in
#' \code{view} passed on to \code{\link[dplyr:filter]{dplyr::filter}()}.
#'
#' @return A mistyR view composition with filtered spatial units from all views.
#'
#' @seealso <\code{\link[dplyr:dplyr_data_masking]{data-masking}}>.
#'
#' @family view manipulation functions
#'
#' @examples
#' # Create a view composition with an intraview and filter
#'
#' library(dplyr)
#'
#' # get the expression data
#' data("synthetic")
#' expr <- synthetic[[1]] %>% select(-c(row, col, type))
#'
#' # compose
#' misty.views <- create_initial_view(expr)
#'
#' # select only the first 10 spatial units and preview
#' misty.views %>%
#'   filter_views(1:10) %>%
#'   str()
#'
#' # select only the units where the expression of ligA is larger than 0.5
#' # and preview
#' misty.views %>%
#'   filter_views(NA, "intraview", ligA > 0.5) %>%
#'   str()
#' @export
filter_views <- function(current.views, rows, view = "intraview", ...) {
  if (!is.numeric(rows)) {
    assertthat::assert_that(length(view) == 1,
      msg = "Please select only a single view for filtering."
    )
    assertthat::assert_that(view %in% names(current.views),
      msg = "The requested view cannot be found in the current view composition."
    )
    assertthat::assert_that(view != "misty.uniqueid",
      msg = "Filtering based on the unique id is not possible."
    )

    # check that filter expressions are logical
    dplyr:::check_filter(dplyr:::dplyr_quosures(...))

    toslice <- misty.views[[view]][["data"]] %>%
      dplyr::transmute(...) %>%
      apply(1, purrr::reduce, `&`) %>%
      which()
  } else {
    toslice <- rows
  }

  purrr::modify_if(
    current.views,
    ~ length(.x) > 1,
    ~ purrr::modify_at(.x, "data", ~ dplyr::slice(.x, toslice))
  )
}


#' Select a subset of markers in a view
#'
#' @inheritParams filter_views
#'
#' @param view the name of the view to select markers for.
#' @param ... one or more \link[dplyr:dplyr_tidy_select]{select} expressions
#' \code{\link[dplyr:select]{dplyr::select}()} for the specified \code{view}.
#'
#' @return A mistyR view composition with selected markers in \code{view}.
#'
#' @seealso <\code{\link[dplyr:dplyr_tidy_select]{tidy-select}}>.
#'
#' @family view manipulation functions
#'
#' @examples
#' # Create a view composition with an intraview and select
#'
#' library(dplyr)
#'
#' # get the expression data
#' data("synthetic")
#' expr <- synthetic[[1]] %>% select(-c(row, col, type))
#'
#' # compose
#' misty.views <- create_initial_view(expr)
#'
#' # select markers from the intraview not starting with lig and preview
#' misty.views %>%
#'   select_markers("intraview", !starts_with("lig")) %>%
#'   str()
#' @export
select_markers <- function(current.views, view = "intraview", ...) {
  assertthat::assert_that(length(view) == 1,
    msg = "Please select only a single view for marker selection."
  )
  assertthat::assert_that(view %in% names(current.views),
    msg = "The requested view cannot be found in the current view composition."
  )
  assertthat::assert_that(view != "misty.uniqueid",
    msg = "Marker selection based on the unique id is not possible."
  )

  purrr::modify_at(
    current.views, view,
    function(x, ...) {
      purrr::modify_at(
        x, "data",
        function(xdata, ...) dplyr::select(xdata, ...),
        ...
      )
    }, ...
  )
}


#' Rename view in a view composition
#'
#' @inheritParams filter_views
#'
#' @param old.name old name of the view.
#' @param new.name new name of the view.
#' @param new.abbrev new abbreviated name.
#'
#' @return A mistyR view composition with a renamed view.
#'
#' @family view manipulation functions
#'
#' @examples
#' view1 <- data.frame(marker1 = rnorm(100, 10, 2), marker2 = rnorm(100, 15, 3))
#' view2 <- data.frame(marker1 = rnorm(100, 10, 5), marker2 = rnorm(100, 15, 5))
#'
#' misty.views <- create_initial_view(view1) %>%
#'   add_views(create_view("originalname", view2, "on"))
#' str(misty.views)
#'
#' # rename and preview
#' misty.views %>%
#'   rename_view("originalname", "renamed", "rn") %>%
#'   str()
#' @export
rename_view <- function(current.views, old.name,
                        new.name, new.abbrev = new.name) {
  assertthat::assert_that(old.name %in% names(current.views),
    msg = "The requested view cannot be found in the current view composition."
  )

  old.view <- purrr::pluck(current.views, old.name)
  old.view$abbrev <- new.abbrev
  changed.view <- list(old.view)
  names(changed.view) <- new.name

  current.views %>%
    rlist::list.remove(old.name) %>%
    append(changed.view)
}
