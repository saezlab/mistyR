library(FNN)
source("multiview.R")

plan(multiprocess, workers=15)

m1 <- read_csv("data/melanoma/33466/33466POST.csv", col_types = cols())

intensities <- m1 %>% select(starts_with("Cell"))
#negative values (-Inf) in Hoechst2, pERK, AXL and MITF
intensities <- as.data.frame(apply(intensities,2,function(i) i <- pmax(0,i)))
positions <- m1 %>% select(CellId, contains("position"))

write_delim(intensities[,-1],"data/melanoma/33466/expressions.txt", delim = " ")
write_csv(positions[,-1],"data/melanoma/33466/positions.txt")

neighborhood <- m1 %>% select(CellId, Percent_Touching, contains("neighbo",ignore.case = T))
stats <- m1 %>% select(CellId, 43:51)



dists <- distances(as.data.frame(positions[,-1]))
l <- (2*400)^2

expr <- intensities %>% select(-CellId)

#time the local multiprocess (4) generation of para view
t0 <- Sys.time()
para.view <- seq(nrow(expr)) %>% future_map_dfr(~data.frame(t(colSums(expr[-.x,] *  exp(-(dists[,.x][-.x]^2)/l)))), .progress = T)
t1 <- Sys.time()

print(t1-t0)

colnames(para.view) <- paste0("p",colnames(para.view))
write_rds(para.view, paste0("para.view.",l,".rds"))


#generate juxta from neighborhood information

t0 <- Sys.time()
juxta.view <- seq(nrow(expr)) %>% future_map_dfr(function(cellid){
  all.neighbors <- neighborhood %>% slice(cellid) %>% select(starts_with("neighbour")) %>% unlist
  all.neighbors <- all.neighbors[all.neighbors > 0]
  
  #maybe
  #The percent_touching neighbors with the smallest disntance to the cell are actually touching
  #percent.touching <- neighborhood %>% slice(cellid) %>% pull(Percent_Touching)
  if(length(all.neighbors) == 0){
    juxta.cell <- (expr %>% slice(1)) * 0 %>% t
  } else {
    juxta.cell <- expr %>% slice(all.neighbors) %>% colSums %>%  t
  }
  
  data.frame(juxta.cell)
}, .progress = TRUE)
t1 <- Sys.time()                   

print(t1-t0)

colnames(juxta.view) <- paste0("j",colnames(juxta.view))
write_rds(juxta.view, "juxta.view.rds")


data <- list.dirs("./data/melanoma", recursive=FALSE)

l <- (2*400)^2
calculateImportances(data,data.l=l,neighbor.thr = NA)

#for comparison
corrplot::corrplot(cor(expr),type = "upper", diag = F)

