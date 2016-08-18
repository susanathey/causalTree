\name{causalForest}
\alias{causalForest}
%\alias{causalTreecallback}
\title{
  Causal Effect Regression and Estimation Forests (Tree Ensembles)
}
\description{
  Build a random forest by fitting a user selected number of \code{causalTree} models to get an ensemble of \code{rpart} objects.
}
\usage{
causalForest(formula, data, treatment,  
                         na.action, 
                         split.Rule, split.Honest, split.Bucket, bucketNum, bucketMax, cv.option, cv.Honest, minsize, 
                         propensity, control, split.alpha, cv.alpha,
                         sample.size.total, sample.size.train.frac,
                         mtry, nodesize, num.trees,
                         cost, weights)
}

    
\arguments{
  \item{formula}{a \link{formula}, with a response and features but no interaction
    terms.  If this a a data frome, that is taken as the model frame
    (see \code{\link{model.frame}).}
  }
  
  \item{data}{an optional data frame that includes the variables
    named in the formula.} 
  
  \item{weights}{optional case weights.}
  
  \item{treatment}{a vector that indicates the treatment status of each observation. 1 represents treated and 0 represents control.  Only binary treatment supported in this version. }

  \item{na.action}{the default action deletes all observations for which
    \code{y} is missing, but keeps those in which one or more predictors
    are missing.}
  
    
    \item{split.Rule}{causalTree splitting options, one of \code{"TOT"}, \code{"CT"}, \code{"fit"}, \code{"tstats"}, four splitting rules in \pkg{causalTree}.  Note that the \code{"tstats"} alternative does not have an associated cross-validation method \code{cv.option}; see Athey and Imbens (2016)
for a discussion.  Note further that \code{split.Rule} and \code{cv.option} can mix and match.} 
    
        
    \item{split.Honest}{boolean option, \code{TRUE} or \code{FALSE}, note that the forest code always calls the honest version of causal tree irrespective of the value of this flag.
}

    \item{split.Bucket}{boolean option, \code{TRUE} or \code{FALSE}, used to specify whether to apply the discrete method in splitting the tree. If set as \code{TRUE}, in splitting a node, the observations in a leaf will be be partitioned into buckets, with each bucket containing \code{bucketNum} treated and \code{bucketNum} control units, and where observations are ordered prior to partitioning. Splitting will take place by bucket.  }
    
    \item{bucketNum}{number of observations in each bucket when set \code{split.Bucket} = \code{TRUE}.  However, the code will override
this choice in order to guarantee that there are at least \code{minsize} and at most \code{bucketMax} buckets.}
    \item{bucketMax}{Option to choose maximum number of buckets to use in splitting when set \code{split.Bucket} = \code{TRUE}, \code{bucketNum} can change by choice of \code{bucketMax}.}
    
    \item{cv.option}{cross validation options, one of \code{"TOT"}, \code{"matching"}, \code{"CT"}, \code{"fit"}, four cross validation methods in \pkg{causalTree}.  There is no \code{cv.option} for the \code{split.Rule} \code{"tstats"}; see Athey and Imbens (2016) for discussion.}
    
    \item{cv.Honest}{boolean option, \code{TRUE} or \code{FALSE}, only used for \code{cv.option} as \code{"CT"} or \code{"fit"}, to specify whether to apply honest risk evalation function in cross validation. If set \code{TRUE}, use honest risk function, otherwise use adaptive risk function in cross validation.  If set \code{FALSE}, the user choice of \code{cv.alpha} will be set to 1.  If set \code{TRUE}, \code{cv.alpha}
will default to 0.5, but the user choice of \code{cv.alpha} will be respected.  Note that honest cv estimates within-leaf variances and may perform better with larger leaf sizes and/or small number of cross-validation sets.}
    
    \item{minsize}{in order to split, each leaf must have at least \code{minsize} treated cases and \code{minsize} control cases. The default value is set as 2.}
    
    \item{propensity}{propensity score used in \code{"TOT"} splitting and \code{"TOT"}, honest \code{"CT"} cross validation methods. The default value is the proportion of treated cases in all observations.  In this implementation, the propensity score is a constant for the whole
dataset.  Unit-specific propensity scores are not supported; however, the user may use inverse propensity scores as case weights if desired.}
    
  \item{control}{a list of options that control details of the
    \code{rpart} algorithm.  See \code{\link{rpart.control}}.}
    
  \item{split.alpha}{scale parameter between 0 and 1, used in splitting risk evaluation function for \code{"CT"}. When \code{split.Honest = FALSE}, \code{split.alpha} will be set as 1.  For \code{split.Rule}=\code{"tstats"}, if \code{split.Honest}=\code{TRUE}, \code{split.alpha} is used in calculating the risk function, which determines the order of pruning in cross-validation.}
  
  \item{cv.alpha}{scale paramter between 0 and 1, used in cross validation risk evaluation function for \code{"CT"} and \code{"fit"}.  When
\code{cv.Honest = FALSE}, \code{cv.alpha} will be set as 1.}
  
  \item{cost}{a vector of non-negative costs, one for each variable in
    the model. Defaults to one for all variables. These are scalings to
    be applied when considering splits, so the improvement on splitting
    on a variable is divided by its cost in deciding which split to
    choose.}
  
  \item{sample.size.total}{Sample size used to build each tree in the forest (sampled randomly with replacement)}
  
  \item{sample.size.train.frac}{Fraction of the sample size used for building each tree (training). For eexample,  if the sample.size.total is 1000 and frac =0.5 then, 500 samples will be used to build the tree and the other 500 samples will be used the evaluate the tree. }
  
  \item{mtry}{Number of data features used to build a tree (This variable is not used presently)}
  
  \item{nodesize}{Minimum number of observations for treated and control cases in one leaf node}
  
  \item{num.trees}{Number of trees to be built in the causal forest}
    
  \item{\dots}{arguments to \code{\link{rpart.control}} may also be
    specified in the call to \code{causalForest}.  They are checked against the
    list of valid arguments.
	The parameter \code{minsize} is implemented differently in \code{causalTree} than in {rpart}; we require a minimum of \code{minsize}
	treated observations and a minimum of \code{minsize} control observations in each leaf.}
}

\details{
CausalForest builds an ensemble of CausalTrees (See Athey and Imbens, \emph{Recursive Partitioning for Heterogeneous Causal
Effects} (2016)), by repeated random sampling of the data with replacement. A causal forest object is a list of trees. To predict, call R's predict function with new test data and the causalForest object (estimated on the training data) obtained after calling the causalForest function. During the prediction phase, the average value over all tree predictions is returned as the final prediction.
 CausalTree differs from \code{rpart} function from \pkg{rpart} package in splitting rules and cross validation methods. Please check Athey and Imbens, \emph{Recursive Partitioning for Heterogeneous Causal
Effects} (2016) and Stefan Wager and Susan Athey, \emph{Estimation and Inference of Heterogeneous Treatment Effects using Random Forests
} for more details.
}

\value{
  An object of class \code{rpart}.  See \code{\link{rpart.object}}.
}

\references{
  Breiman L., Friedman J. H., Olshen R. A., and Stone, C. J. (1984)
  \emph{Classification and Regression Trees.}
  Wadsworth.
  
  Athey, S and G Imbens (2016)  \emph{Recursive Partitioning for Heterogeneous Causal Effects}.  \link{http://arxiv.org/abs/1504.01132}

Wager,S and Athey, S (2015) \emph{Estimation and Inference of Heterogeneous Treatment Effects using Random Forests} \link{http://arxiv.org/abs/1510.04342}
  
}

\seealso{
  \code{\link{causalTree}}
  \code{\link{honest.causalTree}},
  \code{\link{rpart.control}}, \code{\link{rpart.object}},
  \code{\link{summary.rpart}}, \code{\link{rpart.plot}}
}

\examples{
cf <- causalForest(as.formula(paste("y~",f)), data=dataTrain, treatment=dataTrain$w, 
                         split.Rule="CT", split.Honest=T,  split.Bucket=F, bucketNum = 5,
                         bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, 
                        split.alpha = 0.5, cv.alpha = 0.5,
                         sample.size.total = floor(nrow(dataTrain) / 2), sample.size.train.frac = .5,
                         mtry = ceiling(ncol(dataTrain)/3), nodesize = 3, num.trees= 5) 

cfpredtest <- predict(cf, newdata=dataTest, type="vector")
plot(dataTest$tau_true,cfpredtest)

}
\keyword{random forests, tree, causal effects}
