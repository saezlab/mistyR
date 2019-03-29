setwd("/net/data.isilon/ag-saez/bq_jtanevski/spatial")

source("multiview.R")

l <- as.numeric(commandArgs(trailingOnly = T)[1])^2
from <- as.numeric(commandArgs(trailingOnly = T)[2])
to <- as.numeric(commandArgs(trailingOnly = T)[3])

m1 <- read_csv("melanoma/33466POST.csv", col_types = cols())

intensities <- m1 %>% select(starts_with("Cell"))

# negative values (-Inf) in Hoechst2, pERK, AXL and MITF
intensities <- as.data.frame(apply(intensities, 2, function(i) i <- pmax(0, i)))
positions <- m1 %>% select(CellId, contains("position"))

neighborhood <- m1 %>% 
  select(CellId, Percent_Touching, contains("neighbo", ignore.case = T))
stats <- m1 %>% select(CellId, 43:51)

juxta.view <- read_rds("melanoma/juxta.view.rds")

plan(multiprocess, workers = (to - from + 1))

melanoma.views <- create_initial_view(table = intensities[, -1], unique.id = "334466POST") %>%
  add_views(create_view("juxtacrine", "juxta", juxta.view)) %>%
  add_paracrine_view(positions[, -1], l)

estimate_importances(melanoma.views, "melanoma/misty.results/",
                     42, target.subset = seq(from, to))

#alternatively run

plan(list(tweak(multiprocess, workers = (to - from + 1)), 
          tweak(multiprocess, workers = 10)))

estimate_improvement(melanoma.views, paste0("melanoma/misty.results.", l),
                     42, 10, target.subset = seq(from, to)) 