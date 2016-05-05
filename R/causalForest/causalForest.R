library(causalTree)

init.causalForest <- function(x, y, w, num.trees) { 
  num.obs <- nrow(x)
  print("init")
  print(num.obs)
  print(num.trees)
  trees <- vector("list", num.trees) 
  inbag <- matrix(0, num.obs, num.trees) 
  causalForest <- list(trees = trees, x = x, y = y, w = w, ntree = num.trees, inbag = inbag) 
  class(causalForest) <- "causalForest" 
} 

causalForest <- function(X, Y, W, num.trees, sample.size = floor(length(Y) / 10), 
                         mtry = ceiling(ncol(X)/3), nodesize = 1, split.Bucket=F,
                         bucketMax=100, bucketNum=5) {
  
  if (any(is.na(X)) || any(is.na(Y)) || any(is.na(W))) {
    stop("There are missing values in the input.")
  }
  
  num.obs <-nrow(X)
  print(c("num.obs", num.obs))
  causalForest.honest <- init.causalForest(x=X, y=Y, w=W, num.trees=num.trees)
  sample.size <- min(sample.size, floor(num.obs / 2))
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    full.idx <- sample.int(num.obs, 2 * sample.size, replace = FALSE)
    train.idx <- full.idx[1:sample.size]
    reestimation.idx <- full.idx[sample.size + (1:sample.size)]
    
    tree.DF = data.frame(Y = Y, X = X)
    
    tree.honest <- honest.causalTree(Y~., data = tree.DF[train.idx,], 
                                     treatment = W[train.idx], 
                                     est_data=tree.DF[reestimation.idx], est_treatment=W[reestimation.idx],
                                     split.Rule="CT", split.Honest=T, split.Bucket=split.Bucket, 
                                     bucketNum = bucketNum, 
                                     bucketMax = bucketMax, cv.option="CT", cv.Honest=T, 
                                     minsize = nodesize, 
                                     split.alpha = 0.5, cv.alpha = 0.5, xval=0, 
                                     HonestSampleSize=nrow(reestimation.idx), cp=0)
    
    causalForest.honest$trees[[tree.index]] <- tree.honest
    causalForest.honest$inbag[full.idx, tree.index] <- 1
  }
  
  return(causalForest.honest)
}