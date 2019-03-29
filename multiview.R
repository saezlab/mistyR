library(dplyr)
library(purrr)
library(furrr)
library(readr)
library(stringr)
library(tibble)
library(caret)
library(randomForest)
library(deldir)
library(distances)
library(digest)
library(rlist)
library(assertthat)

plan(multiprocess)

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

## Number crunching functions

build_model <- function(views, target, seed = 42, cached = TRUE) {
  set.seed(seed)

  cache.location <- paste0(
    ".misty.temp", .Platform$file.sep,
    views[["misty.uniqueid"]]
  )

  if (!dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  }

  expr <- views[["intracellular"]][["data"]]

  # use the non formula signature of randomForest for memory efficiency
  target.vector <- expr %>% pull(target)

  # returns a list of models
  model.views <- views %>%
    list.remove(c("misty.uniqueid")) %>%
    map(function(view) {
      model.view.cache.file <-
        paste0(
          cache.location, .Platform$file.sep,
          "model.", view[["abbrev"]], ".", target, ".rds"
        )

      if (file.exists(model.view.cache.file) & cached) {
        model.view <- read_rds(model.view.cache.file)
      } else {
        target.index <- match(target, colnames(view[["data"]]))
        model.view <- randomForest(
          x = view[["data"]] %>% select(-target.index),
          y = target.vector, ntree = 100
        )

        if (cached) {
          write_rds(model.view, model.view.cache.file)
        }
      }

      return(model.view)
    })

  # make oob predictions
  oob.predictions <- model.views %>%
    map(~ predict(.x)) %>%
    list.cbind() %>%
    as_tibble() %>%
    add_column(!!target := target.vector)

  # train lm on above
  combined.views <- lm(as.formula(paste0(target, "~.")), oob.predictions)

  # make final.model an object from class misty.model?
  final.model <- list(
    meta.model = combined.views,
    model.views = model.views
  )

  return(final.model)
}


# improvement estimation
estimate_improvement <- function(views, results.folder = "MVResults",
                                 seed = 42, folds = 10, target.subset = NULL) {
  
  if (!dir.exists(results.folder)) dir.create(results.folder, recursive = T)

  expr <- views[["intracellular"]][["data"]]

  header <- "target intra.RMSE intra.R2 multi.RMSE multi.R2"

  write(header, file = paste0(
    results.folder, .Platform$file.sep,
    "performance.txt"
  ))


  targets <- switch(class(target.subset),
    "numeric" = colnames(expr)[target.subset],
    "integer" = colnames(expr)[target.subset],
    "character" = target.subset,
    "NULL" = colnames(expr),
    NULL
  )


  # nested futures! a proper topology should be defined using plan()
  # e.g. plan(list(tweak(multiprocess, workers = 2),
  #               tweak(multiprocess, workers = 4)))
  # for handling 2 targets times 4 folds in parralel
  targets %>% future_map_chr(function(target) {
    set.seed(seed)
    
    target.vector <- expr %>% pull(target)

    test.folds <- createFolds(seq(nrow(expr)), k = folds)

    performance.metrics <- test.folds %>% future_map_dfr(function(fold) {
      # remove test from views
      train.views <- views %>%
        list.remove(c("misty.uniqueid")) %>%
        map(function(view) {
          list(abbrev = view[["abbrev"]], data = view[["data"]] %>%
            slice(-fold))
        }) %>%
        append(list(misty.uniqueid = views[["misty.uniqueid"]]))

      test.views <- views %>%
        list.remove(c("misty.uniqueid")) %>%
        map(function(view) {
          list(abbrev = view[["abbrev"]], data = view[["data"]] %>% slice(fold))
        })

      model.trained <- build_model(train.views, target, seed, FALSE)


      all.predictions <- model.trained[["model.views"]] %>%
        map2(test.views, function(model, view) {
          target.index <- match(target, colnames(view[["data"]]))
          predict(model, view[["data"]] %>%
            select(-target.index) %>%
            mutate(!!target := target.vector[fold]))
        })


      intra.prediction <- all.predictions$intracellular
      multi.view.prediction <- predict(
        model.trained[["meta.model"]],
        as_tibble(all.predictions) %>% mutate(!!target := target.vector[fold])
      )


      intra.RMSE <- RMSE(intra.prediction, target.vector[fold])
      intra.R2 <- R2(intra.prediction, target.vector[fold],
        formula = "traditional"
      )

      multi.RMSE <- RMSE(multi.view.prediction, target.vector[fold])
      multi.R2 <- R2(multi.view.prediction, target.vector[fold],
        formula = "traditional"
      )

      tibble(
        intra.RMSE = intra.RMSE, intra.R2 = intra.R2,
        multi.RMSE = multi.RMSE, multi.R2 = multi.R2
      )
    })

    performance.metrics.summary <- performance.metrics %>% colMeans()

    write(paste(target, paste(performance.metrics.summary, collapse = " ")),
      file = paste0(results.folder, .Platform$file.sep, "performance.txt"),
      append = TRUE
    )

    return(target)
  }, .progress = TRUE)
}


# intracellular should always exist and come first,
# that is why we treat it separately from the other views
estimate_importances <- function(views, results.folder = "MVResults",
                                 seed = 42, target.subset = NULL) {
  if (!dir.exists(results.folder)) dir.create(results.folder, recursive = T)

  view.abbrev <- views %>%
    list.remove(c("misty.uniqueid")) %>%
    map_chr(~ .x[["abbrev"]])


  header <- str_glue("target intercept {views} p.intercept {p.views}",
    views = paste0(view.abbrev, collapse = " "),
    p.views = paste0("p.", view.abbrev, collapse = " "),
    .sep = " "
  )

  expr <- views[["intracellular"]][["data"]]

  write(header, file = paste0(
    results.folder, .Platform$file.sep,
    "coefficients.txt"
  ))

  targets <- switch(class(target.subset),
    "numeric" = colnames(expr)[target.subset],
    "integer" = colnames(expr)[target.subset],
    "character" = target.subset,
    "NULL" = colnames(expr),
    NULL
  )

  targets %>% future_map_chr(function(target) {
    target.model <- build_model(views, target, seed)

    combined.views <- target.model[["meta.model"]]

    model.summary <- summary(combined.views)

    # coefficient values and p-values
    coeff <- c(model.summary$coefficients[, 1], model.summary$coefficients[, 4])

    write(paste(target, paste(coeff, collapse = " ")),
      file = paste0(
        results.folder, .Platform$file.sep,
        "coefficients.txt"
      ), append = TRUE
    )

    # raw importances
    target.model[["model.views"]] %>% walk2(
      view.abbrev,
      function(model.view, abbrev) {
        model.view.imps <- importance(model.view, type = 2)
        imps <- tibble(
          target = rownames(model.view.imps),
          imp = model.view.imps
        )
        write_csv(imps,
          path = paste0(
            results.folder, .Platform$file.sep,
            "RFimportances_", target, "_", abbrev, ".txt"
          )
        )
      }
    )

    return(target)
  }, .progress = TRUE)
}


## Reporting functions
# TODO
