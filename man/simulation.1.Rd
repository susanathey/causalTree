\name{simulation.1}
\alias{simulation.1}
\title{simulation data from Susan and Guido's paper}
\description{
  The \code{simulation.1} data frame has 500 rows and 12 columns, giving 10 x features, corresponding y value, and treatment status for each 500 samples. 
}
\usage{
simulation.1
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{x1},...,\code{x10}}{
      numeric vectors giving the 10 feature values, each follows standard normal distribution.
    }
    \item{\code{y}}{
      numeric vectors giving the response, generated from normal distribution depending on x.
    }
    \item{\code{treatment}}{
    treatment status vector, 1 represents treated and 0 represents controlled.
    }
   }
}
\source{
 One simulation case in Susan and Guido's paper.
}

\examples{
fit <- causalTree(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, 
data = simulation.1, treatment = simulation.1$treatment, split.option = "CT", 
cv.option = "matching", cp = 0, minsize = 5, minbucket = 5)

opcp <- fit$cptable[,1][which.min(fit$cptable[,4])]

opfit <- prune(fit, opcp)

rpart.plot(opfit)
}
\keyword{datasets}
