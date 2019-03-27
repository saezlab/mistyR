library(dplyr)
library(purrr)
library(furrr)
library(readr)
library(stringr)
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
initial_view <- function(table, unique.id = NULL) {
  init.list <- list(intracellular = list(abbrev = "intra", data = table))

  misty.uniqueid <- ifelse(is.null(unique.id), digest(table, "md5"), unique.id)

  view <- append(init.list, list(misty.uniqueid = misty.uniqueid))

  # create cache
  cache.location <- paste0(".misty.temp", .Platform$file.sep, view[["misty.uniqueid"]])

  if (!dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  }

  return(view)
}


create_view <- function(name, abbrev = name, data) {
  new.list <- list(list(abbrev = abbrev, data = data))
  names(new.list)[1] <- name
  return(new.list)
}



add_views <- function(current.views, new.views) {
  # if current views is emtpy fail
  assert_that(length(current.views) >= 1,
    !is.null(current.views[["intracellular"]]),
    !is.null(current.views[["misty.uniqueid"]]),
    msg = "Intracellular view is missing."
  )

  assert_that(is.list(new.views), msg = "The new views are not in a list.")

  assert_that(length(new.views %>% unlist(recursive = F)) %% 2 == 0,
              msg = "The new view is malformed. Consider using create_view().")


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
                msg = "The new view should have the same number  of rows as the intracellular view."
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

  cache.location <- paste0(".misty.temp", .Platform$file.sep, current.views[["misty.uniqueid"]])

  juxta.cache.file <- paste0(cache.location, .Platform$file.sep, "juxta.view.", neighbor.thr, ".rds")

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

  return(current.views %>% add_views(create_view("juxtacrine", paste0("juxta.", neighbor.thr), juxta.view)))
}


add_paracrine_view <- function(current.views, positions, l) {
  dists <- distances(as.data.frame(positions))
  expr <- current.views[["intracellular"]][["data"]]

  cache.location <- paste0(".misty.temp", .Platform$file.sep, current.views[["misty.uniqueid"]])
  para.cache.file <- paste0(cache.location, .Platform$file.sep, "para.view.", l, ".rds")

  if (file.exists(para.cache.file)) {
    cat("Paracrine view retrieved from cache\n")
    para.view <- read_rds(para.cache.file)
  }
  else {
    para.view <- seq(nrow(expr)) %>%
      future_map_dfr(~ data.frame(t(colSums(expr[-.x, ] * exp(-(dists[, .x][-.x]^2) / l)))))
    write_rds(para.view, para.cache.file)
  }

  return(current.views %>% add_views(create_view("paracrine", paste0("para.", l), para.view)))
}


clear_cache <- function(singleid = NULL) {
  if (is.null(singleid)) {
    if (unlink(".misty.temp", recursive = TRUE) == 0) {
      cat("Cache cleared\n")
    } else {
      cat("Failed to clear cache\n")
    }
  } else {
    if (unlink(paste0(".misty.temp", .Platform$file.sep, singleid), recursive = TRUE) == 0) {
      cat("Cache cleared\n")
    } else {
      cat("Failed to clear cache\n")
    }
  }
}


# rewrite below to work with custom views

## Number crunching functions

# improvement estimation
estimate_improvement <- function(views, results.folder = "MVResults",
                                   seed = 42, folds = 10) {

}


# intracellular should always exist and come first, that is why we treat it separately from the other views
estimate_importances <- function(views, results.folder = "MVResults",
                                  seed = 42) {
  set.seed(seed)

  if (!dir.exists(results.folder)) dir.create(results.folder, recursive = T)

  cache.location <- paste0(".misty.temp", .Platform$file.sep, views[["misty.uniqueid"]])


  if (!dir.exists(cache.location)) {
    dir.create(cache.location, recursive = TRUE, showWarnings = TRUE)
  } else {
    cat("Trying to retreive models from cache")
  }

  view.abbrev <- views %>%
    list.remove(c("misty.uniqueid", "intracellular")) %>%
    map_chr(~ list.extract(.x, "abbrev"))

  header <- str_glue("target intercept intra {views} p.intercept p.intra {p.views}",
    views = paste0(view.abbrev, collapse = " "), p.views = paste0("p.", view.abbrev, collapse = " "),
    .sep = " "
  )

  write(header, file = paste0(results.folder, .Platform$file.sep, "coefficients.txt"))

  expr <- views[["intracellular"]][["data"]]

  colnames(expr) %>% future_map(function(target) {

    # use the non formula signature of randomForest for memory efficiency
    target.vector <- expr %>% pull(target)

    model.intra.cache.file <-
      paste0(cache.location, .Platform$file.sep, "model.intra.", target, ".rds")

    if (file.exists(model.intra.cache.file)) {
      model.intra <- read_rds(model.intra.cache.file)
    } else {
      model.intra <- randomForest(
        x = expr %>% select(-target),
        y = target.vector, ntree = 100
      )

      write_rds(model.intra, model.intra.cache.file)
    }

    # returns a list of models
    model.views <- views %>%
      list.remove(c("misty.uniqueid", "intracellular")) %>%
      map(function(view) {
        model.view.cache.file <-
          paste0(cache.location, .Platform$file.sep, "model.", view[["abbrev"]], ".", target, ".rds")

        if (file.exists(model.view.cache.file)) {
          model.view <- read_rds(model.intra.cache.file)
        } else {
          model.view <- randomForest(
            x = view[["data"]] %>% select(-target),
            y = target.vector, ntree = 100
          )
          write_rds(model.view, model.view.cache.file)
        }

        return(model.view)
      })

    # make oob predictions
    oob.predictions <- model.views %>%
      map(~ predict(.x)) %>%
      list.cbind() %>%
      as_tibble() %>%
      add_column(intracellular = predict(model.intra), .before = 1) %>%
      add_column(!!target := target.vector)

    # train lm on above
    combined.views <- lm(as.formula(paste0(target, "~.")), oob.predictions)

    model.summary <- summary(combined.views)

    # coefficient values and p-values
    coeff <- c(model.summary$coefficients[, 1], model.summary$coefficients[, 4])

    write(paste(target, paste(coeff, collapse = " ")),
      file = paste0(results.folder, .Platform$file.sep, "coefficients.txt"), append = TRUE
    )

    # test if there are good rownames
    model.intra.imps <- importance(model.intra, type = 2)
    imps <- tibble(target = rownames(model.intra.imps), imp = model.intra.imps)
    write_csv(imps,
      path = paste0(results.folder, .Platform$file.sep, "RFimportances_", target, "_intra.txt")
    )

    # all other views
    model.views %>% walk2(view.abbrev, function(model.view, abbrev) {
      model.view.imps <- importance(model.view, type = 2)
      imps <- tibble(target = rownames(model.view.imps), imp = model.view.imps)
      write_csv(imps,
        path = paste0(results.folder, .Platform$file.sep, "RFimportances_", target, "_", abbrev, ".txt")
      )
    })

    # TODO: The models are saved, create retreival functions?
    return(list(model.intra, model.views))
  })
}

# TODO: Rewrite
# improvementEstimation <- function(data.paths, results.folder = "MVResults", l = 2^10, neighbor.thr = 1) {
#   set.seed(42)
#
#   data.paths %>% walk(function(d) {
#     if (!dir.exists(paste(d, results.folder, sep = "/"))) dir.create(paste(d, results.folder, sep = "/"))
#
#
#     # data specific / should generalize / pass funciton as argument that takes path at input and returns expr and pos at output
#
#     # fin <- read_csv(paste(d,"random1_position_ALLexpression_real_cells.csv",sep="/"), col_types = cols())
#     # pos <- fin %>% select(c(1,2))
#     # colnames(pos) <- c("x","y")
#     # expr <- fin %>% select(-c(1,2))
#
#     # expr$type <- as.factor(expr$type)
#
#     ###
#
#     expr <- read_delim(paste(d, "expressions.txt", sep = "/"), delim = " ", col_types = cols())
#     colnames(expr) <- make.names(colnames(expr))
#
#     pos <- read_csv(paste(d, "positions.txt", sep = "/"), col_types = cols())
#
#
#     # generate global views
#     # to be used in testing scenario where n cells are queried
#
#     # it HAS to be a data.frame and not a tibble
#     pos <- as.data.frame(pos)
#     dists <- distances(pos)
#
#     if (file.exists(paste0(d, "/", results.folder, "/para.view.", l, ".rds"))) {
#       para.view <- read_rds(paste0(d, "/", results.folder, "/para.view.", l, ".rds"))
#     }
#     else {
#       para.view <- seq(nrow(expr)) %>% future_map_dfr(~ data.frame(t(colSums(expr[-.x, ] * exp(-(dists[, .x][-.x]^2) / l)))))
#       colnames(para.view) <- paste0("p", colnames(para.view))
#       write_rds(para.view, paste0(d, "/", results.folder, "/para.view.", l, ".rds"))
#     }
#
#
#     if (file.exists(paste0(d, "/", results.folder, "/juxta.view.rds"))) {
#       juxta.view <- read_rds(paste0(d, "/", results.folder, "/juxta.view.rds"))
#     }
#     else {
#       delaunay <- deldir(pos)
#
#       juxta.view <- seq(nrow(expr)) %>% future_map_dfr(function(cid) {
#         alln <- getneighbors(delaunay, cid)
#         # first quartile
#         actualn <- alln[which(dists[alln, cid] <= neighbor.thr)]
#
#         data.frame(t(colSums(expr[actualn, ])))
#       })
#
#       juxta.view[is.na(juxta.view)] <- 0
#
#       colnames(juxta.view) <- paste0("j", colnames(juxta.view))
#
#       # hack
#       # colnames(juxta.view)[12:15] <- paste0(colnames(juxta.view)[12:15],"type")
#
#       write_rds(juxta.view, paste0(d, "/", results.folder, "/juxta.view.rds"))
#     }
#
#
#
#     write("target r2.single r2c.single rmse.single r2.multi r2c.multi rmse.multi", file = paste0(d, "/", results.folder, "/improvementEstimation_", l, ".txt"))
#
#     colnames(expr) %>% walk(function(target) {
#
#       # calculate wexpr and juxta using the data from the training data only
#       # test on wexpr and juxta calculated on the whole training set?
#
#       folds <- createFolds(seq(nrow(expr)), k = 10)
#
#
#       multiview.model.performance <- folds %>% future_map_dfr(function(fold) {
#
#         # 3 models for the target (remove pTarget and jTarget)
#
#         # TODO:
#         # no need to retrain model.intra or model.juxta for each l
#         # save the models for each fold in temp folder with '.' prefix and do cleanup(?)
#         # bonus no need to recalculate juxta.view for all l
#
#
#         expr.local <- expr %>% slice(-fold)
#         model.intra <- randomForest(as.formula(paste0(target, "~.")), expr.local, ntree = 100)
#
#
#         # recalculate juxta and para based on training data only
#
#         delaunay.local <- deldir(as.data.frame(pos[-fold, ]))
#
#         juxta.view.local <- seq(nrow(expr.local)) %>% map_dfr(function(cid) {
#           alln <- getneighbors(delaunay.local, cid)
#           # first quartile
#           actualn <- alln[which(dists[alln, cid] <= neighbor.thr)]
#           data.frame(t(colSums(expr[actualn, ])))
#         })
#
#         juxta.view.local[is.na(juxta.view.local)] <- 0
#         colnames(juxta.view.local) <- colnames(juxta.view)
#         j.target <- paste0("j", target)
#         model.juxta <- randomForest(as.formula(paste0(target, "~.")), cbind(juxta.view.local %>% select(-j.target), expr.local %>% select(target)), ntree = 100)
#
#
#         train.ind <- seq(nrow(pos))[-fold]
#         para.view.local <- seq(nrow(expr.local)) %>% map_dfr(~ data.frame(t(colSums(expr.local[-.x, ] * exp(-(dists[train.ind[-.x], train.ind[.x]]^2) / l)))))
#         colnames(para.view.local) <- colnames(para.view)
#         p.target <- paste0("p", target)
#         model.para <- randomForest(as.formula(paste0(target, "~.")), cbind(para.view.local %>% select(-p.target), expr.local %>% select(target)), ntree = 100)
#
#
#         # oob predictions
#         oob.predictions <- data.frame(i = predict(model.intra), j = predict(model.juxta), p = predict(model.para), expr[-fold, ] %>% select(target))
#
#         # train lm on above
#         combined.views <- lm(as.formula(paste0(target, "~.")), oob.predictions)
#
#         # get parameters from lm and make predictions on test data using global para and juxta data
#         single.view.predictions <- predict(model.intra, expr[fold, ])
#         combined.predictions <- predict(combined.views, data.frame(i = predict(model.intra, expr[fold, ]), j = predict(model.juxta, juxta.view[fold, ]), p = predict(model.para, para.view[fold, ]), expr[fold, ] %>% select(target)))
#
#
#
#         r2.single <- R2(pred = single.view.predictions, obs = expr[fold, ] %>% pull(target), formula = "traditional")
#         r2c.single <- R2(pred = single.view.predictions, obs = expr[fold, ] %>% pull(target), formula = "corr")
#         rmse.single <- RMSE(pred = single.view.predictions, obs = expr[fold, ] %>% pull(target))
#
#         r2.multiview <- R2(pred = combined.predictions, obs = expr[fold, ] %>% pull(target), formula = "traditional")
#         r2c.multiview <- R2(pred = combined.predictions, obs = expr[fold, ] %>% pull(target), formula = "corr")
#         rmse.multiview <- RMSE(pred = combined.predictions, obs = expr[fold, ] %>% pull(target))
#
#         c(r2.single, r2c.single, rmse.single, r2.multiview, r2c.multiview, rmse.multiview)
#       })
#
#       print(paste(d, target, paste(rowMeans(multiview.model.performance), collapse = " ")))
#       write(paste(target, paste(rowMeans(multiview.model.performance), collapse = " ")), file = paste0(d, "/", results.folder, "/improvementEstimation_", l, ".txt"), append = TRUE)
#     })
#   })
# }


## Reporting functions

# TODO?
