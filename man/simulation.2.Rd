\name{simulation.2}
\alias{simulation.2}
\title{simulation data from Susan and Guido's paper}
\description{
  The \code{simulation.2} data frame has 500 rows and 12 columns, giving 10 x features, corresponding y value, and treatment status for each 500 samples. 
}
\usage{
simulation.2
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
 One simulation case in Susan and Guido's paper, similar as \code{simulation.1}.
}

\seealso{
\code{\link{simulation.1}}
}

\keyword{datasets}
