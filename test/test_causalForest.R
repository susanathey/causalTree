library(devtools)
#install.packages("rpart", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("rpart.plot", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("reshape2", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("plyr", dependencies=TRUE, repos='http://cran.us.r-project.org')
library(rpart)
library(rpart.plot)
install_github("susanathey/causalTree", ref="forestCode")
# install_github("swager/randomForestCI")

library("causalTree")
library("randomForestCI")


####################################################################################################
# Generate data 
# parameters for data generating
p <- 10 # number of total covariates
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

nameall <- c( name,  "y", "w", "tau_true")

tau_true <- (1-2*w)*(y_ - y)

ntr <- round(.333*n)
nest <- round(.333*n)
ntest <- n - ntr - nest

dataTrain <- data.frame(X[1:ntr,], y[1:ntr], w[1:ntr], tau_true[1:ntr])
dataEst <- data.frame(X[(ntr+1):(ntr+nest),], y[(ntr+1):(ntr+nest)], w[(ntr+1):(ntr+nest)], tau_true[(ntr+1):(ntr+nest)])
dataTest <- data.frame(X[(ntr+nest+1):n,], y[(ntr+nest+1):n], w[(ntr+nest+1):n], tau_true[(ntr+nest+1):n])

names(dataTrain)=nameall
names(dataEst)=nameall
names(dataTest)=nameall



########################## Do analysis


num.trees.temp = ntr

split.Bucket.temp=F
minsize.temp=25
split.Rule.temp = "CT"
cv.option.temp = "CT"
split.Honest.temp = T
cv.Honest.temp = T
split.alpha.temp = .5
cv.alpha.temp = .5

# test out a causal tree just for fun
ct <- causalTree(as.formula(paste("y~",paste(f))), 
         data=dataTrain, treatment=dataTrain$w, 
         split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
         bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
         split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=0, HonestSampleSize=nest, cp=0)
ctpredtest <- predict(ct, newdata=dataTest, type="vector")
plot(dataTest$tau_true,ctpredtest)

# estimate a propensity tree
outcomename = "y"
dataTraintemp <- dataTrain
dataTraintemp$treattreat <- dataTraintemp$w
names(dataTraintemp)[names(dataTraintemp)==outcomename] <- "temptemp"

names(dataTraintemp)[names(dataTraintemp)=="w"] <- outcomename


#one options: estimate the propensity tree with anova so that it will be type "anova" when we re-estimate
#here: replace elements of the rpart object to make it look like anova tree, so that we'll be able to properly predict with it later, etc.
tree.propensity <- rpart(as.formula(paste("y~",f)), data=dataTraintemp, method="class", 
                         control=rpart.control(cp=0, minbucket=minsize.temp*2))

# make it look like a method="anova" tree 
tree.propensity$method <- "anova"
tree.propensity$frame$yval2 <- NULL
tree.propensity$functions$print <- NULL

# switch the names back in the data frame so that when we estimate treatment effects, will have the right outcome variables
names(dataTraintemp)[names(dataTraintemp)==outcomename] <- "w"
names(dataTraintemp)[names(dataTraintemp)=="temptemp"] <- outcomename
pt <- estimate.causalTree(object=tree.propensity,data=dataTraintemp, treatment=dataTraintemp$w)
ptpredtest <- predict(pt, newdata=dataTest, type="vector")
plot(dataTest$tau_true,ptpredtest)


# now estimate a causalForest
cf <- causalForest(as.formula(paste("y~",f)), data=dataTrain, treatment=dataTrain$w, 
                         split.Rule="CT", split.Honest=T,  split.Bucket=F, bucketNum = 5,
                         bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                        split.alpha = 0.5, cv.alpha = 0.5,
                         
                         sample.size.total = floor(nrow(dataTrain) / 2), sample.size.train.frac = .5,
                         mtry = ceiling(ncol(dataTrain)/3), nodesize = 3, num.trees=num.trees.temp) 

cfpredtest <- predict.causalForest(cf, newdata=dataTest, type="vector")
plot(dataTest$tau_true,cfpredtest)

cfpredtrainall <- predict.causalForest(cf, newdata=dataTrain, predict.all = TRUE, type="vector")

# use infJack routine from randomForestCI
cfvar <- infJack(cfpredtrainall$individual, cf$inbag, calibrate = TRUE)
plot(cfvar)

# now estimate a propensityForest
pf <- propensityForest(as.formula(paste("y~",f)), data=dataTrain, treatment=dataTrain$w, 
                   split.Bucket=F, 
                   sample.size.total = floor(nrow(dataTrain) / 2), 
                   mtry = ceiling(ncol(dataTrain)/3), nodesize = 25, num.trees=num.trees.temp) 

pfpredtest <- predict.causalForest(pf, newdata=dataTest, type="vector")
plot(dataTest$tau_true,pfpredtest)

pfpredtrainall <- predict.causalForest(pf, newdata=dataTrain, predict.all = TRUE, type="vector")

pfvar <- infJack(pfpredtrainall$individual, pf$inbag, calibrate = TRUE)
plot(pfvar)

plot(cfpredtest,pfpredtest)
