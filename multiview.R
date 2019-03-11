library(dplyr)
library(purrr)
library(furrr)
library(readr)
library(stringr)
library(caret)
library(randomForest)
library(deldir)
library(distances)

#from a deldir object
getneighbors <- function(ddobj,id){
  union(ddobj$delsgs$ind1[which(ddobj$delsgs$ind2 == id)], ddobj$delsgs$ind2[which(ddobj$delsgs$ind1 == id)])
}

#improvement estimation
improvementEstimation <- function(data.paths, results.folder="MVResults", l = 2^10, neighbor.thr = 1){
  
  set.seed(42)
  
  data.paths %>% walk(function(d){
    
    if(!dir.exists(paste(d, results.folder, sep = "/"))) dir.create(paste(d, results.folder, sep = "/"))
    
    
    #data specific / should generalize / pass funciton as argument that takes path at input and returns expr and pos at output
    
    #fin <- read_csv(paste(d,"random1_position_ALLexpression_real_cells.csv",sep="/"), col_types = cols())
    #pos <- fin %>% select(c(1,2))
    #colnames(pos) <- c("x","y")
    #expr <- fin %>% select(-c(1,2))
    
    #expr$type <- as.factor(expr$type)
    
    ###
    
    expr <- read_delim(paste(d,"expressions.txt",sep="/"), delim = " ", col_types = cols())
    colnames(expr) <- make.names(colnames(expr))
    
    pos <- read_csv(paste(d,"positions.txt",sep="/"), col_types = cols())
    
    
    #generate global views
    #to be used in testing scenario where n cells are queried
    
    #it HAS to be a data.frame and not a tibble
    pos <- as.data.frame(pos)
    dists <- distances(pos)
    
    if(file.exists(paste0(d,"/",results.folder,"/para.view.",l,".rds"))){
      para.view <- read_rds(paste0(d,"/",results.folder,"/para.view.",l,".rds"))
    }
    else{
      para.view <- seq(nrow(expr)) %>% future_map_dfr(~data.frame(t(colSums(expr[-.x,] *  exp(-(dists[,.x][-.x]^2)/l)))))
      colnames(para.view) <- paste0("p",colnames(para.view))
      write_rds(para.view, paste0(d,"/",results.folder,"/para.view.",l,".rds"))
    }
    
    
    if(file.exists(paste0(d,"/",results.folder,"/juxta.view.rds"))){
      juxta.view <- read_rds(paste0(d,"/",results.folder,"/juxta.view.rds"))
    }
    else{
      delaunay <- deldir(pos)
      
      juxta.view <- seq(nrow(expr)) %>% future_map_dfr(function(cid){
        alln <- getneighbors(delaunay,cid)
        #first quartile
        actualn <- alln[which(dists[alln,cid] <= neighbor.thr)]
        
        data.frame(t(colSums(expr[actualn,])))
      })
      
      juxta.view[is.na(juxta.view)] <- 0
      
      colnames(juxta.view) <- paste0("j",colnames(juxta.view))
      
      #hack
      #colnames(juxta.view)[12:15] <- paste0(colnames(juxta.view)[12:15],"type")
      
      write_rds(juxta.view, paste0(d,"/",results.folder,"/juxta.view.rds"))
    }
    
    
    
    write("target r2.single r2c.single rmse.single r2.multi r2c.multi rmse.multi",file = paste0(d, "/", results.folder, "/improvementEstimation_", l, ".txt"))
    
    colnames(expr) %>% walk(function(target){
      
      #calculate wexpr and juxta using the data from the training data only
      #test on wexpr and juxta calculated on the whole training set?
      
      folds <- createFolds(seq(nrow(expr)),k=10)
      
      
      multiview.model.performance <- folds %>% future_map_dfr(function(fold){
        
        #3 models for the target (remove pTarget and jTarget)
        
        #TODO:
        #no need to retrain model.intra or model.juxta for each l
        #save the models for each fold in temp folder with '.' prefix and do cleanup(?)
        #bonus no need to recalculate juxta.view for all l
        
        
        expr.local <- expr %>% slice(-fold)
        model.intra <- randomForest(as.formula(paste0(target,"~.")), expr.local, ntree=100)
        
        
        #recalculate juxta and para based on training data only
        
        delaunay.local <- deldir(as.data.frame(pos[-fold,]))
        
        juxta.view.local <- seq(nrow(expr.local)) %>% map_dfr(function(cid){
          alln <- getneighbors(delaunay.local,cid)
          #first quartile
          actualn <- alln[which(dists[alln,cid] <= neighbor.thr)]
          data.frame(t(colSums(expr[actualn,])))
        })
        
        juxta.view.local[is.na(juxta.view.local)] <- 0
        colnames(juxta.view.local) <- colnames(juxta.view)
        j.target <- paste0("j",target)
        model.juxta <- randomForest(as.formula(paste0(target,"~.")), cbind(juxta.view.local %>% select(-j.target), expr.local %>% select(target)), ntree=100)
        
        
        train.ind <- seq(nrow(pos))[-fold]
        para.view.local <- seq(nrow(expr.local)) %>% map_dfr(~data.frame(t(colSums(expr.local[-.x,] *  exp(-(dists[train.ind[-.x],train.ind[.x]]^2)/l)))))
        colnames(para.view.local) <- colnames(para.view)
        p.target <- paste0("p",target)
        model.para <- randomForest(as.formula(paste0(target,"~.")), cbind(para.view.local %>% select(-p.target), expr.local %>% select(target)), ntree=100)
        
        
        #oob predictions
        oob.predictions <- data.frame(i = predict(model.intra), j = predict(model.juxta), p = predict(model.para), expr[-fold,] %>% select(target))
        
        #train lm on above
        combined.views <- lm(as.formula(paste0(target,"~.")), oob.predictions)
        
        #get parameters from lm and make predictions on test data using global para and juxta data
        single.view.predictions <- predict(model.intra, expr[fold,])
        combined.predictions <- predict(combined.views, data.frame(i = predict(model.intra, expr[fold,]), j = predict(model.juxta, juxta.view[fold,]), p = predict(model.para, para.view[fold,]), expr[fold,] %>% select(target)))
        
        
        
        r2.single <- R2(pred = single.view.predictions, obs = expr[fold, ] %>% pull(target), formula = "traditional")
        r2c.single <- R2(pred = single.view.predictions, obs = expr[fold, ] %>% pull(target), formula = "corr")
        rmse.single <- RMSE(pred = single.view.predictions, obs = expr[fold, ] %>% pull(target))
        
        r2.multiview <- R2(pred = combined.predictions, obs = expr[fold, ] %>% pull(target),formula = "traditional")
        r2c.multiview <- R2(pred = combined.predictions, obs = expr[fold, ] %>% pull(target),formula = "corr")
        rmse.multiview <- RMSE(pred = combined.predictions, obs = expr[fold, ] %>% pull(target))
        
        c(r2.single, r2c.single, rmse.single, r2.multiview, r2c.multiview, rmse.multiview)
        
      })
      
      print(paste(d, target, paste(rowMeans(multiview.model.performance), collapse=" ")))
      write(paste(target, paste(rowMeans(multiview.model.performance), collapse=" ")), file = paste0(d, "/", results.folder, "/improvementEstimation_", l, ".txt"), append=TRUE)
      
    })
    
  })
}


#raw variable importance estimation
calculateImportances <- function(data.paths, results.folder = "MVResults", data.l = 2^10, neighbor.thr = 1){
  
  set.seed(42)
  
  #l should be vector of values, one for each dataset
  if(length(data.l) == 1) data.l <- rep(data.l,length(data.paths))
  stopifnot(length(data.l) == length(data.paths))
  
  #encapsulate in walk
  data.paths %>% iwalk(function(d,ind){
    
    l <- data.l[ind]
    #read or recalculate juxta and para
    expr <- read_delim(paste(d,"expressions.txt",sep="/"), delim = " ",col_types = cols())
    colnames(expr) <- make.names(colnames(expr))
    
    pos <- read_csv(paste(d,"positions.txt",sep="/"), col_names = c("x","y"),col_types = cols())
    
    
    pos <- as.data.frame(pos)
    dists <- distances(pos) 
    
    
    if(file.exists(paste0(d,"/",results.folder,"/para.view.",l,".rds"))){
      para.view <- read_rds(paste0(d,"/",results.folder,"/para.view.",l,".rds"))
    }
    else{
      para.view <- seq(nrow(expr)) %>% future_map_dfr(~data.frame(t(colSums(expr[-.x,] *  exp(-(dists[,.x][-.x]^2)/l)))))
      colnames(para.view) <- paste0("p",colnames(para.view))
      write_rds(para.view, paste0(d,"/",results.folder,"/para.view.",l,".rds"))
    }
    
    
    if(file.exists(paste0(d,"/",results.folder,"/juxta.view.rds"))){
      juxta.view <- read_rds(paste0(d,"/",results.folder,"/juxta.view.rds"))
    }
    else{
      #it HAS to be a data.frame and not a tibble
      delaunay <- deldir(as.data.frame(pos))
      
      juxta.view <- seq(nrow(expr)) %>% future_map_dfr(function(cid){
        alln <- getneighbors(delaunay,cid)
        #first quartile
        actualn <- alln[which(dists[alln,cid] <= neighbor.thr)]
        
        data.frame(t(colSums(expr[actualn,])))
      })
      
      juxta.view[is.na(juxta.view)] <- 0
      
      colnames(juxta.view) <- paste0("j",colnames(juxta.view))
      
      write_rds(juxta.view, paste0(d,"/",results.folder,"/juxta.view.rds"))
    }
    
    write("target intercept intra juxta para p.intercept p.intra p.juxta p.para", file = paste0(d, "/", results.folder, "/coefficients_", l, ".txt"))
    
    colnames(expr) %>% future_map(function(target){
      model.intra <- randomForest(as.formula(paste0(target,"~.")), expr, ntree=100)
      
      j.target <- paste0("j",target)
      model.juxta <- randomForest(as.formula(paste0(target,"~.")), cbind(juxta.view %>% select(-j.target), expr %>% select(target)), ntree=100)
      
      p.target <- paste0("p",target)
      model.para <- randomForest(as.formula(paste0(target,"~.")), cbind(para.view %>% select(-p.target), expr %>% select(target)), ntree=100)
      
      oob.predictions <- data.frame(i = predict(model.intra), j = predict(model.juxta), p = predict(model.para), expr %>% select(target))
      
      #train lm on above
      combined.views <- lm(as.formula(paste0(target,"~.")), oob.predictions)
      
      model.summary <- summary(combined.views)
      
      #coefficient values and p-values
      coeff <- c(model.summary$coefficients[,1], model.summary$coefficients[,4])
      
      imps <- data.frame(target = rownames(importance(model.intra)), imp = importance(model.intra, type=2))
      write_csv(imps, path = paste0(d, "/", results.folder, "/RFimportances_", target, "_", l, "_intra.txt"))
      imps <- data.frame(target = rownames(importance(model.juxta)), imp = importance(model.juxta, type=2))
      write_csv(imps, path = paste0(d, "/", results.folder, "/RFimportances_", target, "_", l, "_juxta.txt"))
      imps <- data.frame(target = rownames(importance(model.para)), imp = importance(model.para))
      write_csv(imps, path = paste0(d, "/", results.folder, "/RFimportances_", target, "_", l, "_para.txt"))
      
      write(paste(target, paste(coeff,collapse = " ")), file = paste0(d, "/", results.folder, "/coefficients_", l, ".txt"), append=TRUE)
    })
    
    print(paste(d,"processed"))
  })
  
}