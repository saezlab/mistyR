## View manipulation functions

# first function to be called
#' Create the basic intrinsic view for MISTy
#'
#' @param table 
#' @param unique.id 
#'
#' @return
#' @export
#'
#' @examples #TBD
create_initial_view <- function(table, unique.id = NULL) {
  init.list <- list(intracellular = list(abbrev = "intra", data = table))

  misty.uniqueid <- ifelse(is.null(unique.id),
    digest::digest(table, "md5"),
    unique.id
  )

  view <- append(init.list, list(misty.uniqueid = misty.uniqueid))

  # create cache
  cache.location <- paste0(
    ".misty.temp", .Platform$file.sep,
    view[["misty.uniqueid"]]
  )

  if (!dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  }

  return(view)
}

# make a misty.view class?
#' Create a view from a data frame
#'
#' @param name 
#' @param data 
#' @param abbrev 
#'
#' @return
#' @export
#'
#' @examples #TBD
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
#' @examples #TBD
add_views <- function(current.views, new.views) {
  assertthat::assert_that(length(current.views) >= 1,
    !is.null(current.views[["intracellular"]]),
    !is.null(current.views[["misty.uniqueid"]]),
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

  view.abbrev <- current.views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map_chr(~ .x[["abbrev"]])

  new.view.abbrev <- new.views %>% purrr::map_chr(~ .x[["abbrev"]])

  assertthat::assert_that(!any(new.view.abbrev %in% view.abbrev),
    msg = "The list of new views contains a duplicate abbreviation."
  )

  new.views %>% walk(function(new.view) {
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
#'
#' @return
#' @export
#'
#' @examples #TBD
add_juxtaview <- function(current.views, positions, neighbor.thr = 15) {
  # from a deldir object
  get_neighbors <- function(ddobj, id) {
    dplyr::union(
      ddobj$delsgs$ind1[which(ddobj$delsgs$ind2 == id)],
      ddobj$delsgs$ind2[which(ddobj$delsgs$ind1 == id)]
    )
  }

  expr <- current.views[["intracellular"]][["data"]]

  cache.location <- paste0(
    ".misty.temp", .Platform$file.sep,
    current.views[["misty.uniqueid"]]
  )

  juxta.cache.file <- paste0(
    cache.location, .Platform$file.sep,
    "juxta.view.", neighbor.thr, ".rds"
  )

  if (file.exists(juxta.cache.file)) {
    cat("Juxtacrine view retrieved from cache\n")
    juxta.view <- readr::read_rds(juxta.cache.file)
  }
  else {
    delaunay <- deldir::deldir(as.data.frame(positions))
    
    juxta.view <- seq(nrow(expr)) %>% furrr::future_map_dfr(function(cid) {
      alln <- get_neighbors(delaunay, cid)
      # suboptimal placement of dists, but makes conflict if out of scope
      # probably due to azy evaluations
      dists <- distances::distances(as.data.frame(positions))
      actualn <- alln[which(dists[alln, cid] <= neighbor.thr)]
      data.frame(t(colSums(expr[actualn, ])))
    })

    readr::write_rds(juxta.view, juxta.cache.file)
  }

  return(current.views %>% add_views(create_view(
    "juxtacrine",
    juxta.view,
    paste0("juxta.", neighbor.thr)
  )))
}


# ncells has priority over nystrom
#' Add paraview to the pipeline
#'
#' @param current.views 
#' @param positions 
#' @param l 
#' @param approx 
#' @param ncells 
#'
#' @return
#' @export
#'
#' @examples #TBD
add_paraview <- function(current.views, positions, l, approx = 1, ncells = NULL) {
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
  expr <- current.views[["intracellular"]][["data"]]

  cache.location <- paste0(
    ".misty.temp", .Platform$file.sep,
    current.views[["misty.uniqueid"]]
  )

  para.cache.file <- paste0(
    cache.location, .Platform$file.sep,
    "para.view.", l, ".rds"
  )

  if (file.exists(para.cache.file)) {
    cat("Paracrine view retrieved from cache\n")
    para.view <- readr::read_rds(para.cache.file)
  }
  else {
    if (is.null(ncells)) {
      if (approx == 1) {
        para.view <- seq(nrow(expr)) %>%
          furrr::future_map_dfr(~ data.frame(t(colSums(expr[-.x, ] *
            exp(-(dists[, .x][-.x]^2) / l)))))
      }
      else {
        if (approx < 1) approx <- base::round(approx * ncol(dists))

        # single Nystrom approximation expert, given RBF with paramter l
        s <- sort(sample.int(n = ncol(dists), size = approx))
        C <- exp(-(dists[, s]^2) / l)

        # pseudo inverse of W
        W.plus <- MASS::ginv(C[s, ])
        # return Nystrom list
        K.approx <- list(s = s, C = C, W.plus = W.plus)

        para.view <- seq(nrow(expr)) %>%
          furrr:future_map_dfr(~ data.frame(t(colSums(expr[-.x, ] * sample_nystrom_row(K.approx, .x)[-.x]))))
      }
    } else {
      para.view <- seq(nrow(expr)) %>%
        furrr::future_map_dfr(function(rowid) {
          knn <- distances::nearest_neighbor_search(dists, ncells + 1, query_indices = rowid)[-1, 1]
          data.frame(t(colSums(expr[knn, ] * exp(-(dists[knn, rowid]^2) / l))))
        })
    }
    readr::write_rds(para.view, para.cache.file)
  }

  return(current.views %>% add_views(create_view(
    paste0("paracrine,", l),
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
#' @examples #TBD
remove_views <- function(current.views, view.names) {
  to.match <- !(view.names %in% c("intracellular", "misty.uniqueid"))
  view.indexes <- match(view.names[to.match], names(current.views))
  current.views %>% rlist::list.remove(view.indexes)
}
