library(devtools)
#install.packages("rpart", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("rpart.plot", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("reshape2", dependencies=TRUE, repos='http://cran.us.r-project.org')
#install.packages("plyr", dependencies=TRUE, repos='http://cran.us.r-project.org')
library(rpart)
library(rpart.plot)
install_github("susanathey/causalTree", ref="modVR1",force=TRUE)
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

ntr <- round(.9*n)
ntest <- n - ntr 

dataTrain <- data.frame(X[1:ntr,], y[1:ntr], w[1:ntr], tau_true[1:ntr])
dataTest <- data.frame(X[(ntr+1):n,], y[(ntr+1):n], w[(ntr+1):n], tau_true[(ntr+1):n])

names(dataTrain)=nameall
names(dataTest)=nameall



########################## Do analysis


num.trees.temp = min(ntr,1000)

split.Bucket.temp=F
minsize.temp=25
split.Rule.temp = "CT"
cv.option.temp = "CT"
split.Honest.temp = T
cv.Honest.temp = T
split.alpha.temp = .5
cv.alpha.temp = .5


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
opcpid <- which.min(tree.propensity$cp[,4])
opcp <- tree.propensity$cp[opcpid,1]
tree.propensity <- prune(tree.propensity, cp = opcp) 

# make it look like a method="anova" tree 
tree.propensity$method <- "anova"
tree.propensity$frame$yval2 <- NULL
tree.propensity$functions$print <- NULL

# switch the names back in the data frame so that when we estimate treatment effects, will have the right outcome variables
names(dataTraintemp)[names(dataTraintemp)==outcomename] <- "w"
names(dataTraintemp)[names(dataTraintemp)=="temptemp"] <- outcomename
pt <- estimate.causalTree(object=tree.propensity,data=dataTraintemp, treatment=dataTraintemp$w)

ptpredtest <- predict(pt, newdata=dataTest, type="vector")
ptpredtrain <- predict(pt, newdata=dataTrain, type="vector")
print(c("mean of ATE treatment effect from propensity tree on Training data", round(mean(ptpredtrain),5)))
plot(dataTest$tau_true,ptpredtest)

ncov_sample<-floor(p/3) #number of covariates (randomly sampled) to use to build tree
ncolx<-p
# now estimate a causalForest
cf <- causalForest(as.formula(paste("y~",f)), data=dataTrain, treatment=dataTrain$w, 
                         split.Rule="CT", split.Honest=T,  split.Bucket=F, bucketNum = 5,
                         bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                        split.alpha = 0.5, cv.alpha = 0.5,
                         sample.size.total = floor(nrow(dataTrain) / 2), sample.size.train.frac = .5,
                         mtry = ceiling(ncol(dataTrain)/3), nodesize = 3, num.trees= 5,ncolx=ncolx,ncov_sample=ncov_sample
                   ) 

cfpredtest <- predict(cf, newdata=dataTest, type="vector")
plot(dataTest$tau_true,cfpredtest)

cfpredtrainall <- predict(cf, newdata=dataTrain, predict.all = TRUE, type="vector")

print(c("mean of ATE treatment effect from causalForest on Training data", round(mean(cfpredtrainall$aggregate),5)))

# use infJack routine from randomForestCI
cfvar <- infJack(cfpredtrainall$individual, cf$inbag, calibrate = TRUE)
plot(cfvar)

# now estimate a propensityForest
pf <- propensityForest(as.formula(paste("y~",f)), data=dataTrain, treatment=dataTrain$w, 
                   split.Bucket=F, 
                   sample.size.total = floor(nrow(dataTrain) / 2), 
                   mtry = ceiling(ncol(dataTrain)/3), nodesize = 25, num.trees=num.trees.temp) 

pfpredtest <- predict(pf, newdata=dataTest, type="vector")
plot(dataTest$tau_true,pfpredtest)

pfpredtrainall <- predict(pf, newdata=dataTrain, predict.all = TRUE, type="vector")
print(c("mean of ATE treatment effect from propensityForest on Training data", round(mean(pfpredtrainall$aggregate),5)))

pfvar <- infJack(pfpredtrainall$individual, pf$inbag, calibrate = TRUE)
plot(pfvar)

plot(cfpredtest,pfpredtest)
