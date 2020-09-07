## View manipulation functions

# first function to be called
#' Create the basic intrinsic view for MISTy
#'
#' @param table
#'
#' @return
#' @export
#'
#' @examples
create_intracellular_view <- function(table) {
  intracellular.view <- list(intracellular = list(abbrev = "intra", data = table))
  return(intracellular.view)
}

get_expression_data <- function(views) {
  return(views[["intracellular"]][["data"]])
}

#' Create a view from a data frame
#'
#' @param name
#' @param data
#' @param abbrev
#'
#' @return
#' @export
#'
#' @examples
create_view <- function(name, data, abbrev = name) {
  new.list <- list(list(abbrev = abbrev, data = data))
  names(new.list)[1] <- name
  return(new.list)
}

#' Add custom views to the current pipeline
#'
#' @param current.views
#' @param new.views
#'
#' @return
#' @export
#'
#' @examples
add_views <- function(current.views, new.views) {
  assertthat::assert_that(length(current.views) >= 1,
    !is.null(current.views[["intracellular"]]),
    msg = "Intracellular view is missing."
  )

  assertthat::assert_that(is.list(new.views),
    msg = "The new views are not in a list or vector."
  )

  assertthat::assert_that(length(new.views %>% unlist(recursive = F)) %% 2 == 0,
    msg = "The new view is malformed. Consider using create_view()."
  )

  assertthat::assert_that(!any(names(new.views) %in% names(current.views)),
    msg = "The list of new views contains a duplicate view name."
  )

  view.abbrev <- current.views %>% purrr::map_chr(~ .x[["abbrev"]])
  new.view.abbrev <- new.views %>% purrr::map_chr(~ .x[["abbrev"]])

  assertthat::assert_that(!any(new.view.abbrev %in% view.abbrev),
    msg = "The list of new views contains a duplicate abbreviation."
  )

  new.views %>% purrr::walk(function(new.view) {
    # check for naming of each element, abbreviation and existance of a table
    assertthat::assert_that(is.list(new.view),
      !is.null(new.view[["abbrev"]]),
      is.character(new.view[["abbrev"]]),
      !is.null(new.view[["data"]]),
      (is.data.frame(new.view[["data"]]) | tibble::is_tibble(new.view[["data"]])),
      msg = "The new view is malformed. Consider using create_view()."
    )

    # check for row compatibility
    assertthat::assert_that(nrow(current.views[["intracellular"]][["data"]]) ==
      nrow(new.view[["data"]]),
    msg = "The new view should have the same number of rows as the intracellular view."
    )
  })

  # update
  return(append(current.views, new.views))
}

#' Add juxtaview to the pipeline
#'
#' @param current.views 
#' @param positions 
#' @param neighbor.thr 
#' @param cached 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
add_juxtaview <- function(current.views, positions, neighbor.thr = 15, cache=NULL, 
                          verbose = TRUE) {
 
  juxta.view.info <- c("juxta.view.", neighbor.thr)
  if(!is.null(cache)) {
    juxta.view <- read_cache_view(cache, juxta.view.info)
  }
  
  if(is.null(cache) || is.null(juxta.view)) {
    if (verbose) message("Computing triangulation")
    delaunay <- deldir::deldir(as.data.frame(positions))
  
    if (verbose) message("Generating juxtaview")
    expression <- get_expression_data(current.views)
    juxta.view <- calculate_juxtaview(expression, positions, delaunay, neighbor.thr, verbose=TRUE)
    
    if(!is.null(cache)) { 
      write_cache_view(cache, juxta.view, juxta.view.info) 
    }
  }

  return(current.views %>% add_views(create_view(
    paste0("juxtaview.", neighbor.thr),
    juxta.view,
    paste0("juxta.", neighbor.thr)
  )))
}

calculate_juxtaview <- function(expression, positions, delaunay, neighbor.thr, verbose=TRUE) {
  
  #ddobj is the result of delaunay triangulation
  get_neighbors <- function(ddobj, id) {
    dplyr::union(
      ddobj$delsgs$ind1[which(ddobj$delsgs$ind2 == id)],
      ddobj$delsgs$ind2[which(ddobj$delsgs$ind1 == id)]
    )
  }
  
  juxta.view <- seq(nrow(expression)) %>% furrr::future_map_dfr(function(cid) {
    alln <- get_neighbors(delaunay, cid)
    # suboptimal placement of dists, but makes conflict if out of scope
    # probably due to lazy evaluations
    dists <- distances::distances(as.data.frame(positions))
    actualn <- alln[which(dists[alln, cid] <= neighbor.thr)]
    data.frame(t(colSums(expression[actualn, ])))
  }, .progress = verbose)
  
  return(juxta.view)
}

# There are two possible approximation methods for paraview: ncells and nystrom. 
# Currently, ncells has priority over nystrom. 
#
#' Add paraview to the pipeline
#'
#' @param current.views 
#' @param positions 
#' @param l 
#' @param approx 
#' @param ncells 
#' @param cached 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
add_paraview <- function(current.views, positions, l, approx = 1, ncells = NULL, cache = NULL, 
                         verbose = TRUE) {
  # K.approx is a list containing C, W.plus and s (indexes of sampled columns)
  sample_nystrom_row <- function(K.approx, k) {

    # transform k into the row index of reordered K.approx
    k.ind <- which(K.approx$s == k)
    if (purrr::is_empty(k.ind)) {
      k.ind <- length(K.approx$s) + k
    }

    cw <- seq(ncol(K.approx$W.plus)) %>%
      purrr::map_dbl(~ K.approx$C[k.ind, ] %*% K.approx$W.plus[, .x])

    cwct <- seq(ncol(t(K.approx$C))) %>%
      purrr::map_dbl(~ cw %*% t(K.approx$C)[, .x])

    # reorder the columns of cwct so that they correspond to the original order
    cwct[c(K.approx$s, seq_along(cwct)[-s])]
  }

  dists <- distances::distances(as.data.frame(positions))
  expression <- get_expression_data(current.views)
  paraview.info <- c("para.view.", l)
  
  if(!is.null(cache)) {
    para.view <- read_cache_view(cache, paraview.info)
  }

  if(is.null(cache) || is.null(para.view)) {
    if (is.null(ncells)) {
      if (approx == 1) {
        if (verbose) message("Generating paraview")
        para.view <- seq(nrow(expression)) %>%
          furrr::future_map_dfr(~ data.frame(t(colSums(expression[-.x, ] *
            exp(-(dists[, .x][-.x]^2) / l)))), 
            .options = furrr::future_options(packages = "distances"))
      }
      else {
        if (approx < 1) approx <- base::round(approx * ncol(dists))

        if (verbose) message("Approximating RBF matrix using the Nystrom method")
        # single Nystrom approximation expert, given RBF with paramter l
        s <- sort(sample.int(n = ncol(dists), size = approx))
        C <- exp(-(dists[, s]^2) / l)

        # pseudo inverse of W
        W.plus <- MASS::ginv(C[s, ])
        # return Nystrom list
        K.approx <- list(s = s, C = C, W.plus = W.plus)

        if (verbose) message("Generating paraview")
        para.view <- seq(nrow(expression)) %>%
          furrr:future_map_dfr(~ data.frame(t(colSums(expression[-.x, ] * sample_nystrom_row(K.approx, .x)[-.x]))))
      }
    } else {
      message("Generating paraview using ", ncells, " nearest neighbors per unit")
      para.view <- seq(nrow(expression)) %>%
        furrr::future_map_dfr(function(rowid) {
          knn <- distances::nearest_neighbor_search(dists, ncells + 1, query_indices = rowid)[-1, 1]
          data.frame(t(colSums(expression[knn, ] * exp(-(dists[knn, rowid]^2) / l))))
        }, .options = furrr::future_options(packages = "distances"))
    }
    
    if (!is.null(cache)) { 
      write_cache_view(cache, para.view, paraview.info) 
    }
  }

  return(current.views %>% add_views(create_view(
    paste0("paraview.", l),
    para.view,
    paste0("para.", l)
  )))
}


#' Remove views from the pipeline
#'
#' @param current.views
#' @param view.names
#'
#' @return
#' @export
#'
#' @examples
remove_views <- function(current.views, view.names) {
  to.match <- !(view.names %in% c("intracellular"))
  view.indexes <- match(view.names[to.match], names(current.views))
  current.views %>% rlist::list.remove(view.indexes)
}
