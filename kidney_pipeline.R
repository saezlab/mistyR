setwd("/net/data.isilon/ag-saez/bq_jtanevski/spatial")

source("multiview.R")
plan(multiprocess, workers = 22)

image <- commandArgs(trailingOnly = T)[1]

d <- readRDS(paste0("kidney/",image,".rds"))

#futures cannot export more than 500MB. Change to 2000 MB.
options(future.globals.maxSize=2000*(1024^2))

# create mcu view based on 100 nearest neigbors
dists <- distances(d[, 24:25])
t <- as.data.frame(table(unique(d$MCU)) - 1)

if(file.exists(paste0("kidney/",image,"_mcu.rds"))){
	mcu.view <- read_rds(paste0("kidney/",image,"_mcu.rds"))
} else {
	# create mcu view based on 100 nearest neigbors
	dists <- distances(d[, 24:25])
	t <- as.data.frame(table(unique(d$MCU)) - 1)
	mcu.view <- seq(nrow(d)) %>%
	future_map(~ table(d$MCU[nearest_neighbor_search(dists, 100, .)[-1, 1]]), .progress = T) %>%
	future_map_dfr(function(cell) {
			mcu.freq <- xtabs(Freq ~ ., rbind(as.data.frame(cell), t))
			to.return <- data.frame(t(as.numeric(mcu.freq)))
			colnames(to.return) <- names(mcu.freq)
			to.return
			}, .progress = T) %>%
	mutate(self.MCU = as.factor(d$MCU))
	write_rds(mcu.view, paste0("kidney/",image,"_mcu.rds"))
}


mpx.views <- create_initial_view(d[, 1:22], unique.id = paste0("4i_",image)) %>%
  add_views(create_view("mcu", data = mcu.view)) %>%
  add_paracrine_view(d[, 24:25], 50^2, ncells = 100) %>%
  add_paracrine_view(d[, 24:25], 100^2, ncells = 100)

estimate_importances(mpx.views, paste0("kidney/results/",image),
                       42,
                      replace = FALSE, sampsize = min(ceiling(.632*nrow(d)), 100000),
                      mtry = max(5, floor(sqrt(ncol(d) - 4))))
