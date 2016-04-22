\name{honest.causalTree}
\alias{honest.causalTree}
%\alias{causalTreecallback}
\title{
  Causal Effect Regression and Estimation Trees: One-step honest estimation
}
\description{
  Fit a \code{causalTree} model to get an honest causal tree, with tree structure built on training sample (including cross-validation) and leaf estimates taken from estimation sample. Return an \code{rpart} object.
}
\usage{
honest.causalTree(formula, data, weights, treatment, subset, 
	est_data, est_weights, est_treatment, est_subset,
	na.action = na.causalTree,
           split.Rule, split.Honest, HonestSampleSize, split.Bucket, bucketNum = 5, 
           bucketMax = 100, cv.option, cv.Honest, minsize = 2L, model = FALSE, x = FALSE,
           y = TRUE, propensity, control, split.alpha = 0.5, cv.alpha = 0.5, cost, \dots)
}

    
\arguments{

  
  \item{est_data}{data frame to be used for leaf estimates; the estimation sample. Must contain the variables used in training the tree.} 
  
  \item{est_weights}{optional case weights for estimation sample}
  
  \item{est_treatment}{treatment vector for estimation sample.  Must be same length as estimation data. A vector indicates the treatment status of the data, 1 represents treated and 0 represents control.  Only binary treatment
			supported in this version. }
  
  \item{est_subset}{optional expression saying that only a subset of the
    rows of the estimation data should be used in the fit of the re-estimated tree.}
  
  \item{...}{All other arguments follow \code{\link{causal.Tree} documentation}}
  
}    
   

\value{
  An object of class \code{rpart}.  See \code{\link{rpart.object}}.
}

\references{
  Breiman L., Friedman J. H., Olshen R. A., and Stone, C. J. (1984)
  \emph{Classification and Regression Trees.}
  Wadsworth.
  
  Athey, S and G Imbens (2016)  \emph{Recursive Partitioning for Heterogeneous Causal Effects}.  \link{http://arxiv.org/abs/1504.01132}
  
}

\seealso{
  \code{\link{causal.Tree} documentation},
  \code{\link{rpart.control}}, \code{\link{rpart.object}},
  \code{\link{summary.rpart}}, \code{\link{print.rpart}}
}

\examples{
tree <- causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                  split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5, 
                  cp = 0, minsize = 20, propensity = 0.5)
                  
opcp <- fit$cptable[,1][which.min(tree$cptable[,4])]

opfit <- prune(tree, opcp)

rpart.plot(opfit)

fittree <- causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                  split.Rule = "fit", cv.option = "fit", split.Honest = T, cv.Honest = T, split.Bucket = T, bucketNum = 5,
		  bucketMax = 200, xval = 10, 
                  cp = 0, minsize = 20, propensity = 0.5)

tstatstree <- causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                  split.Rule = "tstats", cv.option = "CT", cv.Honest = T, split.Bucket = T, bucketNum = 10,
		  bucketMax = 200, xval = 5, 
                  cp = 0, minsize = 20, propensity = 0.5)
}
\keyword{tree, causal effects}
