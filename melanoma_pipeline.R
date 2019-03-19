library(tidyverse)
library(distances)
library(furrr)

m1 <- read_csv("data/melanoma/33466POST.csv", col_types = cols())

intensities <- m1 %>% select(starts_with("Cell"))
positions <- m1 %>% select(CellId, contains("position"))

#The percent_touching neighbors with the smallest disntance to the cell are actually touching
neighborhood <- m1 %>% select(CellId, Percent_Touching, contains("neighbo",ignore.case = T))
stats <- m1 %>% select(CellId, 43:51)

plan(multiprocess, workers=4)

dists <- distances(as.data.frame(positions[,-1]))
l <- (2*400)^2

expr <- intensities %>% select(-CellId)

#time the local multiprocess (4) generation of para view
t0 <- Sys.time()
para.view <- seq(nrow(expr)) %>% future_map_dfr(~data.frame(t(colSums(expr[-.x,] *  exp(-(dists[,.x][-.x]^2)/l)))), .progress = T)
t1 <- Sys.time()

print(t1-t0)

saveRDS(para.view, "para.view.Rds")
                                              
