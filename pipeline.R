#load functions
source("multiview.R")

#configure
data <- list.dirs("./breastcancer", recursive=FALSE)
plan(multiprocess, workers=5)
q25 <- 11.236



#optimize parameter l for each target
#TODO: Can we run (bayesian) optimization here, optimizing the difference in r2? Parameter vector is the value of l for each target
seq(20,400,by=20) %>% walk(~improvementEstimation(data,l=.x^2,neighbor.thr=q25))


#select best l for each target by dataset
results <- data %>% map(function(d){
  seq(20,400,by=20) %>% map(~read_delim(paste0(d,"/MVResults/improvementEstimation_",.x^2,".txt"), skip = 1, col_names = c("target", "r2.single", "r2c.single" , "rmse.single", "r2.multi", "r2c.multi","rmse.multi"), delim = " ", col_types = cols()) %>% transmute(target = target, rmse.imp = (rmse.single - rmse.multi)/rmse.single, r2.imp=pmax(r2.multi,0)-pmax(r2.single,0), r2c.imp=r2c.multi-r2c.single, l = .x^2))
})  

#select based on best improvement over all l
summary <-results %>% map(~reduce(.x, function(accu,x){
  accu$l <- ifelse(x$rmse.imp > accu$rmse.imp, x$l, accu$l)
  accu$r2.imp <- pmax(x$r2.imp, accu$r2.imp)
  accu$rmse.imp <- pmax(x$rmse.imp, accu$rmse.imp)
  accu$r2c.imp <- pmax(x$r2c.imp, accu$r2c.imp)
  accu
}))

#report on improvement in absolute difference in r2

target.imps <- summary %>% map_dfc(~.x[,3]) %>% t
colnames(target.imps) <- summary[[1]]$target

par(mar=c(10,5,2,2)) 

imp.p.value <- apply(target.imps,2,function(x) t.test(x, alternative="g")$p.value)
boxplot(target.imps*100, ylim = c(-8,12), las=2, ylab="Variance explained by cell-cell signaling (%)")

text(x=seq(26),y=rep(11,26),ifelse(imp.p.value<=0.05,"*"," "))

#... and rmse
target.imps <- summary %>% map_dfc(~.x[,2]) %>% t
colnames(target.imps) <- summary[[1]]$target

imp.p.value <- apply(target.imps,2,function(x) t.test(x, alternative="g")$p.value)
boxplot(target.imps*100, las=2, ylab="Improvement in RMSE (%)")

text(x=seq(26),y=rep(26,26),ifelse(imp.p.value<=0.05,"*"," "))


#run importance estimation

#get average l per dataset from the significantly improved targets using RMSE (more liberal)
data.l <- summary %>% map_dbl(~mean(t(.x[imp.p.value<=0.05,5])))
#do it
calculateImportances(data, data.l=data.l, neighbor.thr = q25)




#collect coefficients and calculate percentage of significant per view + distribution of values

#fraction of significant coefficients across datasets
sig.coeff <- data %>% map2(data.l,~read_delim(paste0(.x,"/MVResults/","coefficients_",.y,".txt"), delim= " ", col_types=cols())) %>% reduce(function(accu, coef){
   accu[,2:4] <- accu[,2:4] + as.numeric((coef %>% arrange(target) %>% select(7:9)) <= 0.05)
   accu
}, .init = data.frame(prot=sort(colnames(target.imps)),sigi=0,sigj=0,sigp=0))

sig.coeff <- sig.coeff[order(sig.coeff$prot),]
barplot(as.matrix(t(sig.coeff[imp.p.value[sort(colnames(target.imps))] <= 0.05,2:4]/length(data))), beside = T, names.arg = sig.coeff$prot[imp.p.value[sort(colnames(target.imps))] <= 0.05], las=2, col=c("gray50", "darkseagreen","blue"))


#mean value of the coefficients across datasets
val.coeff <- data %>% map2(data.l,~read_delim(paste0(.x,"/MVResults/","coefficients_",.y,".txt"), delim= " ", col_types=cols())) %>% reduce(function(accu, coef){
  accu[,2:4] <- accu[,2:4] + abs(coef %>% arrange(target) %>% select(3:5)) #* as.numeric((coef %>% arrange(target) %>% select(7:9)) <= 0.05))
  accu
}, .init = data.frame(prot= sort(colnames(target.imps)) ,ci=0, cj=0, cp=0))

val.coeff <- val.coeff[order(val.coeff$prot),]
barplot(as.matrix(t(val.coeff[imp.p.value[sort(colnames(target.imps))] <= 0.05,2:4]/length(data))), beside = F, names.arg = val.coeff$prot[imp.p.value[sort(colnames(target.imps))] <= 0.05], las=2, col=c("gray50", "darkseagreen","blue"))


# valc <- data %>% map(~read_delim(paste0(.x,"/MLResults/coeff_8192.txt"),delim= " ",col_names = F, col_types = cols())) %>% reduce(function(accu, coef){
#   accu[,2:4] <- accu[,2:4] + coef[,3:5] * as.numeric(coef[,7:9]<=0.05)
#   accu
# }, .init = data.frame(prot=colnames(expr),ci=0,cj=0,cp=0))
# 
# #mean
# valc <- valc[order(valc$prot),]
# valc[,2:4] <- valc[,2:4]/sigc[,2:4]
# 
# barplot(as.matrix(t(abs(valc[pvals$p<=0.05,2:4]))),beside = F,names.arg = valc$prot[pvals$p<=0.05], las=2, col=c("gray50", "darkseagreen","blue"))



#for the significantly improved targets, get the intra para and juxta importances
#WARNING: recycling target.imps

intra.importances <- matrix(data = 0, ncol(target.imps), ncol(target.imps))
juxta.importances <- matrix(data = 0, ncol(target.imps), ncol(target.imps))
para.importances <- matrix(data = 0, ncol(target.imps), ncol(target.imps))


data %>% walk2(data.l, function(data.path,l){
  coefficients <- read_delim(paste0(data.path,"/MVResults/","coefficients_",l,".txt"), delim= " ", col_types=cols())
  
  intra.impmap <- matrix(data = NA, ncol(target.imps), ncol(target.imps))
  juxta.impmap <- matrix(data = NA, ncol(target.imps), ncol(target.imps))
  para.impmap <- matrix(data = NA, ncol(target.imps), ncol(target.imps))
  
  colnames(intra.impmap) <- colnames(target.imps)
  rownames(intra.impmap) <- colnames(target.imps)
  intra.impmap <- as.data.frame(intra.impmap)
  
  colnames(juxta.impmap) <- paste0("j",colnames(target.imps))
  rownames(juxta.impmap) <- colnames(target.imps)
  juxta.impmap <- as.data.frame(juxta.impmap)
  
  colnames(para.impmap) <- paste0("p",colnames(target.imps))
  rownames(para.impmap) <- colnames(target.imps)
  para.impmap <- as.data.frame(para.impmap)
  
  #collect importances from results folder, standardize per row and multiply by 1 - p.value from coefficients
  colnames(target.imps) %>% walk(function(prot){
    current.imp <- read_csv(paste0(data.path,"/MVResults/","RFimportances_",prot,"_",l,"_intra.txt"), col_types=cols())
    intra.impmap[prot, current.imp$target] <<- (current.imp$IncNodePurity - mean(current.imp$IncNodePurity))/sd(current.imp$IncNodePurity) * (1 - coefficients[which(coefficients$target == prot),]$p.intra)
    current.imp <- read_csv(paste0(data.path,"/MVResults/","RFimportances_",prot,"_",l,"_juxta.txt"), col_types=cols())
    juxta.impmap[prot, current.imp$target] <<- (current.imp$IncNodePurity - mean(current.imp$IncNodePurity))/sd(current.imp$IncNodePurity) * (1 - coefficients[which(coefficients$target == prot),]$p.juxta)
    current.imp <- read_csv(paste0(data.path,"/MVResults/","RFimportances_",prot,"_",l,"_para.txt"), col_types=cols())
    para.impmap[prot, current.imp$target] <<- (current.imp$IncNodePurity - mean(current.imp$IncNodePurity))/sd(current.imp$IncNodePurity) * (1 - coefficients[which(coefficients$target == prot),]$p.para)
  })
  
  #aggregate into global importances by summing
  intra.importances <<- intra.importances + intra.impmap
  juxta.importances <<- juxta.importances + juxta.impmap
  para.importances <<- para.importances + para.impmap
  
})

#filter the rows for the significantly improved targets
intra.importances %<>% slice(which(imp.p.value <= 0.05))
rownames(intra.importances) <- names(imp.p.value)[which(imp.p.value <= 0.05)]

juxta.importances %<>% slice(which(imp.p.value <= 0.05))
rownames(juxta.importances) <- names(imp.p.value)[which(imp.p.value <= 0.05)]

para.importances %<>% slice(which(imp.p.value <= 0.05))
rownames(para.importances) <- names(imp.p.value)[which(imp.p.value <= 0.05)]


#AWS with p-value weights?


#heatmaps and functions for network extraction/visualization
library(pheatmap)

pheatmap(intra.importances, cluster_rows = F, cluster_cols = F)
pheatmap(juxta.importances, cluster_rows = F, cluster_cols = F)
pheatmap(para.importances, cluster_rows = F, cluster_cols = F)






#Synthetic data

# toComply <- function(d){
#   datain <- read_csv(paste0(d,"/random1_position_ALLexpression_real_cells.csv"))
#   pos <- datain[,1:2]
#   colnames(pos) <- c("x","y")
#   expr <- datain[,-c(1,2)]
#   write_csv(pos,paste0(d,"/positions.txt"))
#   write_delim(expr[,-ncol(expr)], paste0(d,"/expressions.txt"), delim = " ")
# }
# 
# data <- list.dirs("./synthetic", recursive = FALSE)
# plan(multiprocess, workers=5)
# q25 <- 1
# 
# data %>% walk(~toComply(.x))
# 
# c(2,5,seq(10,80,by=10)) %>% walk(~improvementEstimation(data,l=.x^2,neighbor.thr=q25))
