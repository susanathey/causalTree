# This program tests out the features of the causalTree package

#install.packages("devtools")
library(devtools)
#install.packages("rpart", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("rpart.plot", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("reshape2", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("plyr", dependencies=TRUE, repos='http://cran.us.r-project.org')
library(rpart)
library(rpart.plot)
# install_github("susanathey/causalTree",force=TRUE)
library(causalTree)
library(reshape2)
library(plyr)


# Generate data 
# parameters for data generating
p <- 20 # number of total covariates
pt <- 4 # number of covariates affecting treatment effects
py <- 4 # number of covariates affecting outcomes but not treatment effects
asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
n <- 3000 # total size of the dataset
propens <- .5 #treatment probability
sig = .01
treatsize <- .5 # treatment effect size
levsize <- 1

# draw W
w <- rbinom(n, 1, propens)
w_ <- 1-w

# draw X
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# generate treatment effects as function of X
if (p<pt+py) print("error: p>=pt+py required")
tau <- 0
for (iii in 1:pt) {
  tau <- tau + treatsize*pmax(X[,iii],array(0,n))*(2*X[,iii])
}

# generate average value of outcomes
mu <- treatsize*rowSums(X[,1:pt])+levsize*rowSums(X[,(pt+1):(pt+py)])

# generate outcomes as function of treatment status, mu, tau, and noise
y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)
y_ <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)

# create formulas for estimation
# if modifying code for another dataset, need to simply name outcome variable y, treatment variable w, and
# create a formula such as f with the list of x variables separated by "+"
f <- ""
nextx <- ""
if (p>1) {
  for (ii in 1:(p-1)) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
    f <- paste(f, nextx, "+", sep="")
  }
  f <- paste(f, "x", ii+1, sep="")
} else if (p==1) {
  f <- "x1"
}

for (ii in 1:p) {
  nextx <- paste("x",ii, sep="")
  if (ii==1) {name <- nextx}
  if (ii>1) {name <- c(name, nextx)}
}

name <- c( name,  "y", "w", "tau_true")

tau_true <- (1-2*w)*(y_ - y)

ntr <- round(.333*n)
nest <- round(.333*n)
ntest <- n - ntr - nest

# set global parameters
minsize.temp = 25
split.Bucket.temp = T
bucketNum.temp = 5
bucketMax.temp = 100

X<-data.frame(X)
for (tmp1 in 1:ncol(X)){
  xtmp<-X[,tmp1]
  unxtmp<-unique(xtmp)
  if(length(unxtmp)<bucketMax.temp) #convert to factor
    X[,tmp1]<-factor(X[,tmp1])
}

dataTrain <- data.frame(X[1:ntr,], y[1:ntr], w[1:ntr], tau_true[1:ntr])
dataEst <- data.frame(X[(ntr+1):(ntr+nest),], y[(ntr+1):(ntr+nest)], w[(ntr+1):(ntr+nest)], tau_true[(ntr+1):(ntr+nest)])
dataTest <- data.frame(X[(ntr+nest+1):n,], y[(ntr+nest+1):n], w[(ntr+nest+1):n], tau_true[(ntr+nest+1):n])

names(dataTrain)=name
names(dataEst)=name
names(dataTest)=name

tree_honest_prune_list <- vector(mode="list", length=4)
tree_dishonest_prune_list <- vector(mode="list", length=4)

# preselect cross-validation groups to remove randomness in comparing methods
xvalvec = sample(5, nrow(dataTrain), replace=TRUE)
#xvalvec=5


# Do causal tree estimation
split.Rule.temp = "CT" #CT
cv.option.temp = "CT" #CT
split.Honest.temp = T
cv.Honest.temp = T
split.alpha.temp = .5
cv.alpha.temp = .5



#This function is a wrapper for honest causal tree
tree <- honest.causalTree(as.formula(paste("y~",paste(f))),
                  data=dataTrain, treatment=dataTrain$w,
                  est_data=dataEst, est_treatment=dataEst$w,
                  split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp,
                  bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp,
                  split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
#You can still prune as usual; the cptable is the one from training the tree
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp)

# save the results
tree_honest_CT <- tree
tree_honest_CT_prune <- tree_prune

#manually convert honest to dishonest: honest causaltree re-estimation using training data itself
tree <- honest.causalTree(as.formula(paste("y~",paste(f))),
                          data=dataTrain, treatment=dataTrain$w,
                          est_data=dataTrain, est_treatment=dataTrain$w,
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp,
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp,
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)


# get the dishonest version--estimated leaf effects on training sample
tree <- causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                   split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                   bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                   split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 
tree_dishonest_CT_prune <- tree_prune

# we can get the honest tree manually, by estimating the leaf effects on a new sample
tree_honest_CT_prune2 <- estimate.causalTree(object=tree_dishonest_CT_prune,data=dataEst, treatment=dataEst$w)

print(tree_honest_CT_prune)
print(tree_honest_CT_prune2)

#add to list
tree_dishonest_prune_list[[1]] <- tree_dishonest_CT_prune
tree_honest_prune_list[[1]] <- tree_honest_CT_prune



# Other honest splitting rules, cv

#tstats
split.Rule.temp = "tstats"
cv.option.temp = "CT"
split.Honest.temp = F
cv.Honest.temp = T
split.alpha.temp = 1

tree <- honest.causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                          est_data=dataEst, est_treatment=dataEst$w,
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
# prune
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_honest_tstats <- tree
tree_honest_tstats_prune <- tree_prune
tree_honest_prune_list[[2]] <- tree_honest_tstats_prune

#honest fit
split.Rule.temp = "fit"
cv.option.temp = "fit"
split.Honest.temp = T
cv.Honest.temp = T
split.alpha.temp = .5

tree <- honest.causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                          est_data=dataEst, est_treatment=dataEst$w,
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
#You can still prune as usual; the cptable is the one from training the tree
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_honest_fit <- tree
tree_honest_fit_prune <- tree_prune

tree_honest_prune_list[[3]] <- tree_honest_fit_prune

#honest TOT
split.Rule.temp = "TOT"
cv.option.temp = "TOT"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = 1

tree <- honest.causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                          est_data=dataEst, est_treatment=dataEst$w,
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
#You can still prune as usual; the cptable is the one from training the tree
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_honest_TOT <- tree
tree_honest_TOT_prune <- tree_prune
tree_honest_prune_list[[4]] <- tree_honest_TOT_prune


##################################
#Dishonest approaches


# CT
split.Rule.temp = "CT"
cv.option.temp = "CT"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = 1

tree <- causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                   split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                   bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                   split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)

opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_dishonest_CT <- tree
tree_dishonest_CT_prune <- tree_prune
tree_dishonest_prune_list[[1]] <- tree_dishonest_CT_prune

#tstats
split.Rule.temp = "tstats"
cv.option.temp = "CT"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = 1

tree <- causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                   split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                   bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                   split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)

opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_dishonest_tstats <- tree
tree_dishonest_tstats_prune <- tree_prune

tree_dishonest_prune_list[[2]] <- tree_dishonest_tstats_prune

#fit
split.Rule.temp = "fit"
cv.option.temp = "fit"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = .5

tree <- causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                          split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                          bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                          split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)

opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_dishonest_fit <- tree
tree_dishonest_fit_prune <- tree_prune
tree_dishonest_prune_list[[3]] <- tree_dishonest_fit_prune

#TOT
split.Rule.temp = "TOT"
cv.option.temp = "TOT"
split.Honest.temp = F
cv.Honest.temp = F
split.alpha.temp = 1

tree <- causalTree(as.formula(paste("y~",paste(f))), 
                          data=dataTrain, treatment=dataTrain$w, 
                   split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                   bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                   split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalvec, HonestSampleSize=nest, cp=0)
#You can still prune as usual; the cptable is the one from training the tree
opcpid <- which.min(tree$cp[,4])
opcp <- tree$cp[opcpid,1]
tree_prune <- prune(tree, cp = opcp) 

# save the results
tree_dishonest_TOT <- tree
tree_dishonest_TOT_prune <- tree_prune
tree_dishonest_prune_list[[4]] <- tree_dishonest_TOT_prune

# performance on test set
# if modifying for a real dataset, comment out the parts with infeasible MSE since we won't know that
MSEtau_infeas_honest <- array(0,4)
MSEtau_infeas_dishonest <- array(0,4)
MSEtau_TOT_honest <- array(0,4)
MSEtau_TOT_dishonest <- array(0,4)
MSEy_honest <- array(0,4)
MSEy_dishonest <- array(0,4)

ystar = dataTest$y*(dataTest$w-propens)/(propens*(1-propens))
for (i in 1:4) {
  # use the predictions from honest trees on the Test set
  predicthonest = predict(tree_honest_prune_list[[i]],newdata=dataTest,type="vector")
  # use the predictions from dishonest trees on Test set
  predictdishonest = predict(tree_dishonest_prune_list[[i]],newdata=dataTest,type="vector")
  # this is simulated data, so we know the true tau for each individual.  This is infeasible MSE
  MSEtau_infeas_honest[[i]] = mean((predicthonest - dataTest$tau_true)^2, na.rm=T)
  MSEtau_infeas_dishonest[[i]] = mean((predictdishonest - dataTest$tau_true)^2, na.rm=T)
  # in a typical dataset, can use TOT loss function
  MSEtau_TOT_honest[[i]] = mean((predicthonest-ystar)^2, na.rm=T)
  MSEtau_TOT_dishonest[[i]] = mean((predictdishonest-ystar)^2, na.rm=T)
  
  # look at how the different partitions do predicting outcomes within a leaf
  dataTest$leaff <- as.factor(round(predicthonest,4))
  
  #honest estimation of predicted y
  dataEst$leaff <- as.factor(round(predict(tree_honest_prune_list[[i]], newdata=dataEst,type="vector"),4))
  yPredHonestTable <- melt(tapply(dataEst$y, list(dataEst$leaff, dataEst$w), mean), varnames=c("leaff","w"))
  yPredHonestTable <- rename(yPredHonestTable,replace=c("value"="ypredhon"))
  dataTest <- merge(dataTest,yPredHonestTable, by.x=c("leaff", "w"), by.y=c("leaff","w"))
  MSEy_honest[[i]] = mean((dataTest$ypredhon-dataTest$y)^2, na.rm=T)
  
  #dishonest estimation of predicted y
  dataTest$leaffd <- as.factor(round(predictdishonest,4))
  dataTrain$leaffd <- as.factor(round(predict(tree_dishonest_prune_list[[i]], newdata=dataTrain,type="vector"),4))
  yPredDishonestTable <- melt(tapply(dataTrain$y, list(dataTrain$leaffd, dataTrain$w), mean), varnames=c("leaffd","w"))
  yPredDishonestTable <- rename(yPredDishonestTable,replace=c("value"="ypreddishon"))
  dataTest <- merge(dataTest,yPredDishonestTable, by.x=c("leaffd", "w"), by.y=c("leaffd","w"))
  MSEy_dishonest[[i]] = mean((dataTest$ypreddishon-dataTest$y)^2, na.rm=T)
  
  dataTest <- dataTest[, !(names(dataTest) %in% c("ypredhon","ypreddishon"))]
}
  
  
print("infeasible MSE(tau)'s for honest CT, tstats, fit, TOT")
print(MSEtau_infeas_honest)
print("infeasible MSE(tau)'s for dishonest CT, tstats, fit, TOT")
print(MSEtau_infeas_dishonest)
print("TOT MSE(tau)'s for honest CT, tstats, fit, TOT")
print(MSEtau_TOT_honest)
print("TOT MSE(tau)'s for dishonest CT, tstats, fit, TOT")
print(MSEtau_TOT_dishonest)
print("MSE(y)'s for honest CT, tstats, fit, TOT")
print(MSEy_honest)
print("MSE(y)'s for dishonest CT, tstats, fit, TOT")
print(MSEy_dishonest)


# easy trick to get coefficients and standard errors

dataTrain$leaves <- predict(tree_dishonest_CT_prune, newdata=dataTrain, type = 'vector')
dataEst$leaves <- predict(tree_dishonest_CT_prune, newdata=dataEst, type = 'vector')
dataTest$leaves <- predict(tree_dishonest_CT_prune, newdata=dataTest, type = 'vector')


dataTrain$leavesf <- factor(round(dataTrain$leaves,4))
dataEst$leavesf <- factor(round(dataEst$leaves,4))
dataTest$leavesf <- factor(round(dataTest$leaves,4))

# run regressions with indicators for the leaves interacted with the treatment indicator
if (length(levels(dataTrain$leavesf)) == 1){
  
  modelTrain <- lm(y~w, data=dataTrain)
  modelEst <- lm(y~w, data=dataEst)
  modelTest <- lm(y~w, data=dataTest)
  
  summary(modelTrain)
  summary(modelEst)
  summary(modelTest)
  
} else{
  
  modelTrain <- lm(y~-1+leavesf+leavesf*w-w, data=dataTrain)
  modelEst <- lm(y~-1+leavesf+leavesf*w-w, data=dataEst)
  modelTest <- lm(y~-1+leavesf+leavesf*w-w, data=dataTest)
  
  print("Leaf names match estimated treatment effects on training set")
  print(summary(modelTrain))
  print("Estimated treatment effects on estimation set typically more moderate than training set")
  print(summary(modelEst))
  print("Estimated treatment effects on test set typically more moderate than training set")
  print(summary(modelTest))
  
  
  # extract the coefficient vectors which are the leaf treatment effects
  coefnumh <- length(coef(modelEst))
  coefnuml <- length(coef(modelEst))/2 + 1
  
  Train.coeftr <- coef(modelTrain)[coefnuml:coefnumh]
  Est.coeftr <- coef(modelEst)[coefnuml:coefnumh]
  Test.coeftr <- coef(modelTest)[coefnuml:coefnumh]
  
  # calculate leaf probabilities
 
  leafprobEst <- tapply(dataEst$y,list(dataEst$leavesf),length)
  leafprobTrain <- tapply(dataTrain$y,list(dataTrain$leavesf),length)
  leafprobTest <- tapply(dataTest$y,list(dataTest$leavesf),length)
  leafprob <- (leafprobEst + leafprobTrain + leafprobTest)/(nrow(dataEst) + nrow(dataTrain) + 
                                                              nrow(dataTest))
  
  #calculate variance of estimated treatment effects--typically this is higher in the training set, since there is overfitting there
  Train.coefvar <- sum(leafprob * Train.coeftr^2)-(sum(leafprob*Train.coeftr)^2)
  Est.coefvar <- sum(leafprob * Est.coeftr^2)-(sum(leafprob*Est.coeftr)^2)
  Test.coefvar <- sum(leafprob * Test.coeftr^2)-(sum(leafprob*Test.coeftr)^2)
  
  print("Variance of estimated treatment effects: Train, Estimation, Test Sets")
  print("Typically train has higher var--more extreme estimates--due to overfitting")
  print(c(Train.coefvar,Est.coefvar,Test.coefvar))
}      

