init.causalForest <- function(formula, data, treatment, weights=F, cost=F, num.trees,ncov_sample) { 
  num.obs <- nrow(data)
  trees <- vector("list", num.trees)
  inbag <- matrix(0, num.obs, num.trees) 
  cov_sample <- matrix(0,num.trees,ncov_sample)
  inbag.Est <- matrix(0, num.obs, num.trees)
  nameall_sample <- matrix(0,num.trees,ncov_sample+2)#2 end cols for y,w,no tau_true
  fsample<-vector("list",num.trees)
  causalForestobj <- list(trees = trees, formula=formula, data=data, treatment=treatment, weights=weights, cost=cost, ntree = num.trees, inbag = inbag,cov_sample=cov_sample, fsample=fsample,nameall_sample=nameall_sample,inbag.Est=inbag.Est) 
  class(causalForestobj) <- "causalForest" 
  return(causalForestobj)
} 

predict.causalForest <- function(forest, newdata, predict.all = FALSE, type="vector") {
  if (!inherits(forest, "causalForest")) stop("Not a legitimate \"causalForest\" object")  

  individual <- sapply(forest$trees, function(tree.fit) {
    predict(tree.fit, newdata=newdata, type="vector")
  })
  
  #replace sapply with a loop if needed
  print(dim(individual))
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
                         cost=F, weights=F,ncolx,ncov_sample) {
  
  # do not implement subset option of causalTree, that is inherited from rpart but have not implemented it here yet

  num.obs <-nrow(data)
  causalForest.hon <- init.causalForest(formula=formula, data=data, treatment=treatment, weights=weights, cost=cost, num.trees=num.trees,ncov_sample=ncov_sample)
  sample.size <- min(sample.size.total, num.obs)
  train.size <- round(sample.size.train.frac*sample.size)
  est.size <- sample.size - train.size
  

  treatmentdf <- data.frame(treatment)
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)
    train.idx <- full.idx[1:train.size]
    reestimation.idx <- full.idx[(train.size+1):sample.size]
    
    #randomize over the covariates for splitting (both train and reestimation)
    cov_sample<-sample.int(ncolx)
    cov_sample<-cov_sample[1:ncov_sample]
    #modify the y=f(x) equation accordingly for this tree
  fsample<-""
  nextx<-""
  if (ncov_sample>1) {
  for (ii in 1:(ncov_sample-1)) {
  nextx <- paste("x",cov_sample[ii], sep="")
  if (ii==1) {name <- nextx}
  if (ii>1) {name <- c(name, nextx)}
  fsample <- paste(fsample, nextx, "+", sep="")
  }
  fsample <- paste(fsample, "x", cov_sample[ii+1], sep="")
  } else if (ncov_sample==1) {
  fsample <- paste("x",cov_sample[1], sep="")
  }
  #modify the colnames
    
    nameall_sample<-c()
  for (ii in 1:ncov_sample) {
  nextx <- paste("x",cov_sample[ii], sep="")
  if (ii==1) {name <- nextx}
  if (ii>1) {name <- c(name, nextx)}
  }
  nameall_sample <- c( name,  "y", "w") #, "tau_true")

    
    #store this var subset for each tree (need it during testing/predict stage)
    causalForest.hon$cov_sample[tree.index,]<-cov_sample
    #also store the formula & colnames of X for each tree (need it during testing/predict stage)
    causalForest.hon$nameall_sample[tree.index,]<-nameall_sample
    causalForest.hon$fsample[[tree.index]]<-fsample
    
    dataTree <- data.frame(data[train.idx,])
    dataEstim <- data.frame(data[reestimation.idx,])
    

    #pick relevant covariates for tree
    dataTree <- dataTree[,c(cov_sample,(ncolx+1):ncol(dataTree))]
    dataEstim <- dataEstim[,c(cov_sample,(ncolx+1):ncol(dataEstim))]
    
    #change colnames to reflect the sampled cols
    names(dataTree)=nameall_sample
    names(dataEstim)=nameall_sample
    
    
    
    
    #save rdata for debug here, if needed
    formula<-paste("y~",fsample,sep="")
    tree.honest <- honest.causalTree(formula, data = dataTree, 
                                     treatment = treatmentdf[train.idx,], 
                                     est_data=dataEstim, est_treatment=treatmentdf[reestimation.idx,],
                                     split.Rule="CT", split.Honest=T, split.Bucket=split.Bucket, 
                                     bucketNum = bucketNum, 
                                     bucketMax = bucketMax, cv.option="CT", cv.Honest=T, 
                                     minsize = nodesize, 
                                     split.alpha = 0.5, cv.alpha = 0.5, xval=0, 
                                     HonestSampleSize=est.size, cp=0)


    causalForest.hon$trees[[tree.index]] <- tree.honest
    causalForest.hon$inbag[full.idx, tree.index] <- 1
    causalForest.hon$inbag.Est[reestimation.idx, tree.index] <- 1
  }
  
  return(causalForest.hon)
}


propensityForest <- function(formula, data, treatment,  
                         na.action = na.causalTree, 
                         split.Rule="CT", split.Honest=T, split.Bucket=F, bucketNum = 5,
                         bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                         propensity=mean(treatment), control, split.alpha = 0.5, cv.alpha = 0.5,  
                         
                         sample.size.total = floor(nrow(data) / 10), sample.size.train.frac = 1,
                         mtry = ceiling(ncol(data)/3), nodesize = 1, num.trees=nrow(data),ncolx=ncolx,ncov_sample=ncov_sample) {
  
  # do not implement subset option of causalTree, inherited from rpart
  # do not implement weights and costs yet
  
  if(sample.size.train.frac != 1) {
    print("warning: for propensity Forest, sample.size.train.frac should be 1; resetting to 1")
    sample.size.train.frac <- 1
  }
  
  num.obs <-nrow(data)
  

  causalForest.hon <- init.causalForest(formula=formula, data=data, treatment=treatment, num.trees=num.trees, weights=F, cost=F,ncov_sample=ncov_sample)
  sample.size <- min(sample.size.total, num.obs)
  train.size <- round(sample.size.train.frac*sample.size)
  
  treatmentdf <- data.frame(treatment)
  outcomename = as.character(formula[2])
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)
    train.idx <- full.idx[1:train.size]
    
    cov_sample<-sample.int(ncolx)
    cov_sample<-cov_sample[1:ncov_sample]
    #modify the y=f(x) equation accordingly for this tree
    fsample<-""
    nextx<-""
    if (ncov_sample>1) {
      for (ii in 1:(ncov_sample-1)) {
        nextx <- paste("x",cov_sample[ii], sep="")
        if (ii==1) {name <- nextx}
        if (ii>1) {name <- c(name, nextx)}
        fsample <- paste(fsample, nextx, "+", sep="")
      }
      fsample <- paste(fsample, "x", cov_sample[ii+1], sep="")
    } else if (ncov_sample==1) {
      fsample <- "x1"
    }
    #modify the colnames
    
    nameall_sample<-c()
    for (ii in 1:ncov_sample) {
      nextx <- paste("x",cov_sample[ii], sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
    }
    
    nameall_sample <- c( name,"temptemp","y", "tau_true","treattreat")
    nameall_sample_save <- c( name,  "y", "w", "tau_true")
    
     #store this var subset for each tree (need it during testing/predict stage)
    causalForest.hon$cov_sample[tree.index,]<-cov_sample
    #also store the formula & colnames of X for each tree (need it during testing/predict stage)
    causalForest.hon$nameall_sample[tree.index,]<-nameall_sample_save
    causalForest.hon$fsample[[tree.index]]<-fsample
      
    # rename variables as a way to trick rpart into building the tree with all the object attributes considering the outcome variable as named
    # by the input formula, even though the tree itself is trained on w.  Note that we aren't saving out this propensity tree anyway, but if
    # we decided later to try to save out the propensity trees and do something directly with the propensity scores, we would need to do something
    # more tedious like estimate the propensity tree with the original names, and then edit the attributes to replace the treatment variable name
    # with the outcome variable name for the estimate part
    dataTree <- data.frame(data[train.idx,])
    dataTree$treattreat <- treatmentdf[train.idx,]
    names(dataTree)[names(dataTree)==outcomename] <- "temptemp"

    names(dataTree)[names(dataTree)=="treattreat"] <- outcomename
    
    
    
    #sample covariates
    #pick relevant covariates for tree
    dataTree <- dataTree[,c(cov_sample,(ncolx+1):ncol(dataTree))]
    # dataEstim <- dataEstim[,c(cov_sample,(ncolx+1):ncol(dataEstim))]
    
    #change colnames to reflect the sampled cols
    names(dataTree)=nameall_sample
    # names(dataEstim)=nameall_sample
    formula<-paste("y~",fsample,sep="")
    #one options: estimate the propensity tree with anova so that it will be type "anova" when we re-estimate
    #here: replace elements of the rpart object to make it look like anova tree, so that we'll be able to properly predict with it later, etc.
    tree.propensity <- rpart(formula=formula, data=dataTree, method="class", 
                             control=rpart.control(cp=0, minbucket=nodesize))
    
    # make it look like a method="anova" tree 
    tree.propensity$method <- "anova"
    tree.propensity$frame$yval2 <- NULL
    tree.propensity$functions$print <- NULL
    
    # switch the names back in the data frame so that when we estimate treatment effects, will have the right outcome variables
    names(dataTree)[names(dataTree)==outcomename] <- "treattreat"
    names(dataTree)[names(dataTree)=="temptemp"] <- outcomename
    tree.treatment <- estimate.causalTree(object=tree.propensity,data=dataTree, treatment=dataTree$treattreat)
    
    causalForest.hon$trees[[tree.index]] <- tree.treatment
    causalForest.hon$inbag[full.idx, tree.index] <- 1
  }
  
  return(causalForest.hon)
}