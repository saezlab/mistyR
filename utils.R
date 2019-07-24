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
