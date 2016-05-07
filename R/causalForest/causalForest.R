library(causalTree)

init.causalForest <- function(formula, data, treatment, weights=F, cost=F, num.trees) { 
  num.obs <- nrow(data)
  formform <- as.formula(paste(as.character(formula[2]),as.character(formula[1]),as.character(formula[3])))
  print(formform)
  trees <- vector("list", num.trees)
  print(trees)
  inbag <- matrix(0, num.obs, num.trees) 
  print(inbag)
  causalForestobj <- list(trees = trees, formula=formform, data=data, treatment=treatment, weights=weights, cost=cost, ntree = num.trees, inbag = inbag) 
  class(causalForestobj) <- "causalForest" 
  return(causalForestobj)
} 

predict.causalForest <- function(forest, newdata, predict.all = FALSE, type="vector") {
  if (!inherits(forest, "causalForest")) stop("Not a legitimate \"causalForest\" object")  

  individual <- sapply(forest$trees, function(tree.fit) {
    predict(tree.fit, newdata=newdata, type="vector")
  })
  
  aggregate <- rowMeans(individual)
  if (predict.all) {
    list(aggregate = aggregate, individual = individual)
  } else {
    aggregate
  }
}

causalForest <- function(formula, data, treatment,  
                         na.action = na.causalTree, 
                         split.Rule="CT", split.Honest=T, split.Bucket=F, bucketNum = 5,
                         bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                         propensity, control, split.alpha = 0.5, cv.alpha = 0.5,  
                         
                         sample.size.total = floor(nrow(data) / 10), sample.size.train.frac = .5,
                         mtry = ceiling(ncol(data)/3), nodesize = 1, num.trees=nrow(data),
                         cost=F, weights=F) {
  
  # do not implement subset option of causalTree, inherited from rpart

  num.obs <-nrow(data)
  formula2 <- as.formula(paste(as.character(formula[2]),as.character(formula[1]),as.character(formula[3])))
  causalForest.hon <- init.causalForest(formula=formula2, data=data, treatment=treatment, weights=weights, cost=cost, num.trees=num.trees)
  sample.size <- min(sample.size.total, num.obs)
  train.size <- round(sample.size.train.frac*sample.size)
  est.size <- sample.size - train.size
  
  formform <- as.formula(paste(as.character(formula[2]),as.character(formula[1]),as.character(formula[3])))
  treatmentdf <- data.frame(treatment)
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)
    train.idx <- full.idx[1:train.size]
    reestimation.idx <- full.idx[(train.size+1):sample.size]
    
    dataTree <- data.frame(data[train.idx,])
    dataEstim <- data.frame(data[reestimation.idx,])

    tree.honest <- honest.causalTree(formula=formform, data = dataTree, 
                                     treatment = treatmentdf[train.idx,], 
                                     est_data=dataEstim, est_treatment=treatmentdf[reestimation.idx,],
                                     split.Rule="CT", split.Honest=T, split.Bucket=split.Bucket, 
                                     bucketNum = bucketNum, 
                                     bucketMax = bucketMax, cv.option="CT", cv.Honest=T, 
                                     minsize = nodesize, 
                                     split.alpha = 0.5, cv.alpha = 0.5, xval=0, 
                                     HonestSampleSize=est.size, cp=0)

    print(tree.honest)
    causalForest.hon$trees[[tree.index]] <- tree.honest
    causalForest.hon$inbag[full.idx, tree.index] <- 1
  }
  
  return(causalForest.hon)
}


propensityForest <- function(formula, data, treatment,  
                         na.action = na.causalTree, 
                         split.Rule="CT", split.Honest=T, split.Bucket=F, bucketNum = 5,
                         bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                         propensity, control, split.alpha = 0.5, cv.alpha = 0.5,  
                         
                         sample.size.total = floor(nrow(data) / 10), sample.size.train.frac = 1,
                         mtry = ceiling(ncol(data)/3), nodesize = 1, num.trees=nrow(data),
                         cost=F, weights=F) {
  
  # do not implement subset option of causalTree, inherited from rpart
  
  if(sample.size.train.frac != 1) {
    print("warning: for propensity Forest, sample.size.train.frac should be 1; resetting to 1")
    sample.size.train.frac <- 1
  }
  
  num.obs <-nrow(data)
  formula2 <- as.formula(paste(as.character(formula[2]),as.character(formula[1]),as.character(formula[3])))
  causalForest.hon <- init.causalForest(formula=formula2, data=data, treatment=treatment, weights=weights, cost=cost, num.trees=num.trees)
  sample.size <- min(sample.size.total, num.obs)
  train.size <- round(sample.size.train.frac*sample.size)
  
  formform <- as.formula(paste(as.character(formula[2]),as.character(formula[1]),as.character(formula[3])))
  treatmentdf <- data.frame(treatment)
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)
    train.idx <- full.idx[1:train.size]
    
    dataTree <- data.frame(data[train.idx,])
    dataEstim <- data.frame(data[reestimation.idx,])
    
    tree.propensity <- rpart()
    
    tree.honest <- honest.causalTree(formula=formform, data = dataTree, 
                                     treatment = treatmentdf[train.idx,], 
                                     est_data=dataEstim, est_treatment=treatmentdf[reestimation.idx,],
                                     split.Rule="CT", split.Honest=T, split.Bucket=split.Bucket, 
                                     bucketNum = bucketNum, 
                                     bucketMax = bucketMax, cv.option="CT", cv.Honest=T, 
                                     minsize = nodesize, 
                                     split.alpha = 0.5, cv.alpha = 0.5, xval=0, 
                                     HonestSampleSize=est.size, cp=0)
    
    print(tree.honest)
    causalForest.hon$trees[[tree.index]] <- tree.honest
    causalForest.hon$inbag[full.idx, tree.index] <- 1
  }
  
  return(causalForest.hon)
}