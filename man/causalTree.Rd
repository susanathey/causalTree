\name{causalTree}
\alias{causalTree}
%\alias{causalTreecallback}
\title{
  Causal Effect Regression and Estimation Trees
}
\description{
  Fit a \code{causalTree} model to get an \code{rpart} object
}
\usage{
causalTree(formula, data, weights, treatment, subset, na.action = na.causalTree,
           split.Rule, split.Honest, HonestSampleSize, split.Bucket, bucketNum = 5, 
           bucketMax = 100, cv.option, cv.Honest, minsize = 2L, model = FALSE, x = FALSE,
           y = TRUE, propensity, control, split.alpha = 0.5, cv.alpha = 0.5, cost, \dots)
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
  
  \item{subset}{optional expression saying that only a subset of the
    rows of the data should be used in the fit.}
  
  \item{na.action}{the default action deletes all observations for which
    \code{y} is missing, but keeps those in which one or more predictors
    are missing.}
  
    
    \item{split.Rule}{causalTree splitting options, one of \code{"TOT"}, \code{"CT"}, \code{"fit"}, \code{"tstats"}, four splitting rules in \pkg{causalTree}.  Note that the \code{"tstats"} alternative does not have an associated cross-validation method \code{cv.option}; see Athey and Imbens (2016)
for a discussion.  Note further that \code{split.Rule} and \code{cv.option} can mix and match.} 
    
        
    \item{split.Honest}{boolean option, \code{TRUE} or \code{FALSE}, used for \code{split.Rule} as \code{"CT"} or \code{"fit"}. If set as \code{TRUE}, do honest splitting, with default \code{split.alpha} = 0.5; if set as \code{FALSE}, do adaptive splitting with \code{split.alpha} = 1.  The user choice of \code{split.alpha} will be ignored if \code{split.Honest} is set as \code{FALSE}, but will be respected
if set to \code{TRUE}.  For \code{split.Rule}=\code{"TOT"}, there is no honest splitting option and
the parameter \code{split.alpha} does not matter.  For \code{split.Rule}=\code{"tstats"}, a value of \code{TRUE} enables use of \code{split.alpha} in calculating the risk function, which determines the order of pruning in cross-validation. Note also that causalTree function
returns the estimates from the training data, no matter what the value of \code{split.Honest} is; the tree must be re-estimated to get the honest estimates using \code{estimate.causalTree}. The wrapper function \code{honest.CausalTree}
does honest estimation in one step and returns a tree.}
    
    \item{HonestSampleSize}{number of observations anticipated to be used in honest re-estimation after building the tree. This enters the
risk function used in both splitting and cross-validation.}
    
    \item{split.Bucket}{boolean option, \code{TRUE} or \code{FALSE}, used to specify whether to apply the discrete method in splitting the tree. If set as \code{TRUE}, in splitting a node, the observations in a leaf will be be partitioned into buckets, with each bucket containing \code{bucketNum} treated and \code{bucketNum} control units, and where observations are ordered prior to partitioning. Splitting will take place by bucket.  }
    
    \item{bucketNum}{number of observations in each bucket when set \code{split.Bucket} = \code{TRUE}.  However, the code will override
this choice in order to guarantee that there are at least \code{minsize} and at most \code{bucketMax} buckets.}
    \item{bucketMax}{Option to choose maximum number of buckets to use in splitting when set \code{split.Bucket} = \code{TRUE}, \code{bucketNum} can change by choice of \code{bucketMax}.}
    
    \item{cv.option}{cross validation options, one of \code{"TOT"}, \code{"matching"}, \code{"CT"}, \code{"fit"}, four cross validation methods in \pkg{causalTree}.  There is no \code{cv.option} for the \code{split.Rule} \code{"tstats"}; see Athey and Imbens (2016) for discussion.}
    
    \item{cv.Honest}{boolean option, \code{TRUE} or \code{FALSE}, only used for \code{cv.option} as \code{"CT"} or \code{"fit"}, to specify whether to apply honest risk evalation function in cross validation. If set \code{TRUE}, use honest risk function, otherwise use adaptive risk function in cross validation.  If set \code{FALSE}, the user choice of \code{cv.alpha} will be set to 1.  If set \code{TRUE}, \code{cv.alpha}
will default to 0.5, but the user choice of \code{cv.alpha} will be respected.  Note that honest cv estimates within-leaf variances and may perform better with larger leaf sizes and/or small number of cross-validation sets.}
    
    \item{minsize}{in order to split, each leaf must have at least \code{minsize} treated cases and \code{minsize} control cases. The default value is set as 2.}
  
  \item{model}{if logical: keep a copy of the model frame in the result?
    If the input value for \code{model} is a model frame (likely from an
    earlier call to the \code{causalTree} function), then this frame is used
    rather than constructing new data.}
  \item{x}{keep a copy of the \code{x} matrix in the result.}
  
  \item{y}{keep a copy of the dependent variable in the result.  If
    missing and \code{model} is supplied this defaults to \code{FALSE}.}
    
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
  
  \item{\dots}{arguments to \code{\link{rpart.control}} may also be
    specified in the call to \code{causalTree}.  They are checked against the
    list of valid arguments.  An example of a commonly set parameter would be \code{xval}, which sets the number of cross-validation samples.
	The parameter \code{minsize} is implemented differently in \code{causalTree} than in {rpart}; we require a minimum of \code{minsize}
	treated observations and a minimum of \code{minsize} control observations in each leaf.}
}

\details{
 CausalTree differs from \code{rpart} function from \pkg{rpart} package in splitting rules and cross validation methods. Please check Athey and Imbens, \emph{Recursive Partitioning for Heterogeneous Causal
Effects} (2016) for more details.
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
  \code{\link{honest.causalTree}},
  \code{\link{rpart.control}}, \code{\link{rpart.object}},
  \code{\link{summary.rpart}}, \code{\link{rpart.plot}}
}

\examples{
tree <- causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                  split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5, 
                  cp = 0, minsize = 20, propensity = 0.5)
                  
opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]

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
