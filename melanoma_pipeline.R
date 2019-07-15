setwd("/net/data.isilon/ag-saez/bq_jtanevski/spatial")

source("multiview.R")

data <- commandArgs(trailingOnly = T)[1]
from <- as.numeric(commandArgs(trailingOnly = T)[2])
to <- as.numeric(commandArgs(trailingOnly = T)[3])

m1 <- read_csv(paste0("melanoma/", data, ".csv"), col_types = cols())

# select measured proteins and remove Hoechst
intensities <- m1 %>% select(starts_with("Cell"), -contains("Hoechst"))

# deal with negative values (-Inf)
intensities <- as.data.frame(apply(intensities, 2, function(i) i <- pmax(0, i)))

positions <- m1 %>% select(CellId, contains("position"))

neighborhood <- m1 %>%
  select(CellId, Percent_Touching, contains("neighbo", ignore.case = T))
stats <- m1 %>% select(CellId, 43:51)

rm(m1)

#futures cannot export more than 500MB. Change to 2000 MB.
options(future.globals.maxSize=2000*(1024^2))

# Start parallel processing
plan(multiprocess, workers = (to - from + 1))

if (file.exists(paste0("melanoma/", data, ".juxta.view.rds"))) {
  juxta.view <- read_rds(paste0("melanoma/", data, ".juxta.view.rds"))
} else {
  expr <- intensities %>% select(-CellId)
  # generate juxta.view from neighborhood
  juxta.view <- seq(nrow(expr)) %>% future_map_dfr(function(cellid) {
    all.neighbors <- neighborhood %>%
      slice(cellid) %>%
      select(starts_with("neighbour")) %>%
      unlist()
    all.neighbors <- all.neighbors[all.neighbors > 0]

    if (length(all.neighbors) == 0) {
      juxta.cell <- (expr %>% slice(1)) * 0 %>% t()
    } else {
      juxta.cell <- expr %>%
        slice(all.neighbors) %>%
        colSums() %>%
        t()
    }

    data.frame(juxta.cell)
  }, .progress = TRUE)

  write_rds(juxta.view, paste0("melanoma/", data, ".juxta.view.rds"))
}

# create views

melanoma.views <- create_initial_view(table = intensities[, -1], unique.id = data) %>%
  add_views(create_view("juxtacrine", "juxta", juxta.view)) %>%
  add_paracrine_view(positions[, -1], l=400^2, ncells=1000) %>%
  add_paracrine_view(positions[, -1], l=800^2, ncells=1000)

 estimate_importances(melanoma.views, paste0("melanoma/triplet_results/",data),
                     42, target.subset = seq(from, to), 
                    replace = FALSE, sampsize = min(ceiling(.632*nrow(intensities)), 100000), 
                    mtry = max(5, floor(sqrt(ncol(intensities) - 1))))

# alternatively run

# plan(list(tweak(multiprocess, workers = (to - from + 1)),
#          tweak(multiprocess, workers = 10)))

# estimate_improvement(melanoma.views, paste0("melanoma/misty.results.", l),
#                     42, 10, target.subset = seq(from, to))
