## View manipulation functions

# first function to be called
create_initial_view <- function(table, unique.id = NULL) {
  init.list <- list(intracellular = list(abbrev = "intra", data = table))
  
  misty.uniqueid <- ifelse(is.null(unique.id), digest(table, "md5"), unique.id)
  
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
create_view <- function(name, abbrev = name, data) {
  new.list <- list(list(abbrev = abbrev, data = data))
  names(new.list)[1] <- name
  return(new.list)
}



add_views <- function(current.views, new.views) {
  assert_that(length(current.views) >= 1,
              !is.null(current.views[["intracellular"]]),
              !is.null(current.views[["misty.uniqueid"]]),
              msg = "Intracellular view is missing."
  )
  
  assert_that(is.list(new.views),
              msg = "The new views are not in a list or vector."
  )
  
  assert_that(length(new.views %>% unlist(recursive = F)) %% 2 == 0,
              msg = "The new view is malformed. Consider using create_view()."
  )
  
  assert_that(!any(names(new.views) %in% names(current.views)),
              msg = "The list of new views contains a duplicate view name."
  )
  
  view.abbrev <- current.views %>%
    list.remove(c("misty.uniqueid")) %>%
    map_chr(~ .x[["abbrev"]])
  
  new.view.abbrev <- new.views %>% map_chr(~ .x[["abbrev"]])
  
  assert_that(!any(new.view.abbrev %in% view.abbrev),
              msg = "The list of new views contains a duplicate abbreviation."
  )
  
  new.views %>% walk(function(new.view) {
    # check for naming of each element, abbreviation and existance of a table
    assert_that(is.list(new.view),
                !is.null(new.view[["abbrev"]]),
                is.character(new.view[["abbrev"]]),
                !is.null(new.view[["data"]]),
                (is.data.frame(new.view[["data"]]) | is_tibble(new.view[["data"]])),
                msg = "The new view is malformed. Consider using create_view()."
    )
    
    # check for row compatibility
    assert_that(nrow(current.views[["intracellular"]][["data"]]) == nrow(new.view[["data"]]),
                msg = "The new view should have the same number of rows as the intracellular view."
    )
  })
  
  # update
  return(append(current.views, new.views))
}

add_juxtacrine_view <- function(current.views, positions, neighbor.thr = 15) {
  # from a deldir object
  get_neighbors <- function(ddobj, id) {
    union(
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
    juxta.view <- read_rds(juxta.cache.file)
  }
  else {
    delaunay <- deldir(as.data.frame(positions))
    dists <- distances(as.data.frame(positions))
    
    juxta.view <- seq(nrow(expr)) %>% future_map_dfr(function(cid) {
      alln <- get_neighbors(delaunay, cid)
      # first quartile
      actualn <- alln[which(dists[alln, cid] <= neighbor.thr)]
      
      data.frame(t(colSums(expr[actualn, ])))
    })
    
    write_rds(juxta.view, juxta.cache.file)
  }
  
  return(current.views %>% add_views(create_view(
    "juxtacrine",
    paste0("juxta.", neighbor.thr),
    juxta.view
  )))
}


add_paracrine_view <- function(current.views, positions, l) {
  dists <- distances(as.data.frame(positions))
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
    para.view <- read_rds(para.cache.file)
  }
  else {
    para.view <- seq(nrow(expr)) %>%
      future_map_dfr(~ data.frame(t(colSums(expr[-.x, ] *
                                              exp(-(dists[, .x][-.x]^2) / l)))))
    write_rds(para.view, para.cache.file)
  }
  
  return(current.views %>% add_views(create_view(
    "paracrine",
    paste0("para.", l),
    para.view
  )))
}


remove_views <- function(current.views, view.names) {
  to.match <- !(view.names %in% c("intracellular", "misty.uniqueid"))
  view.indexes <- match(view.names[to.match], names(current.views))
  current.views %>% list.remove(view.indexes)
}