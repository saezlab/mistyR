create_cache <- function(data, cache.id=NULL){
  
  cache.id <- ifelse(is.null(cache.id),
                    digest::digest(data, "md5"),
                    cache.id
  )
  
  # create cache
  cache.location <- paste0(
    ".misty.temp", .Platform$file.sep,
    cache.id
  )
  
  if (!dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  } else {
    message(paste0("Note: Cache directory with id", cache.id, " already exists."))
  }
  
  cache = list(cach.id=cache.id, location=cache.location)
  return(cache)
}

read_cache_view <- function(cache, view.name) {
  view.name <- paste0(view.name, collapse="")
  cache.file <- paste0(
    cache$location, .Platform$file.sep,
    paste0(c(view.name, ".rds"), collapse="")
  )
  
  print(cache.file)
  if (file.exists(cache.file)) {
    message(paste(view.name, "view retrieved from cache\n"))
    view <- readr::read_rds(cache.file)
    return(view)
  } 
  
  return(NULL)
}

write_cache_view <- function(cache, view.data, view.name) { 
  cache.file <- paste0(
    cache$location, .Platform$file.sep,
    paste0(c(view.name, ".rds"), collapse="")
  )
  
  readr::write_rds(view.data, cache.file)
  return(cache.file)
}

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