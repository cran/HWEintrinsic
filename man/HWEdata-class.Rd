\name{HWEdata-class}
\Rdversion{1.1}
\docType{class}
\alias{HWEdata-class}
\alias{HWEdata}
\alias{initialize,HWEdata-method}

\title{Class "HWEdata". Data specification for the Hardy-Weinberg Testing Problem Using the Intrinsic Prior Approach.}
\description{
  This class encapsulates the data specification for a Bayesian objective analysis via intrinsic priors of the Hardy-Weinberg Testing Problem as described in Consonni et al. (2011).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HWEdata", data)}, where \code{data} are the data in vector form.
}
\section{Slots}{
  \describe{
    \item{\code{data.mat}:}{Object of class \code{"matrix"}; data in matrix form.}
    \item{\code{size}:}{Object of class \code{"numeric"}; number of alleles included in the data.}
    \item{\code{data.vec}:}{Object of class \code{"numeric"}; data in vector form.}
  }
}
\section{Methods}{
  \describe{
    %\item{initialize}{\code{signature(.Object = "HWEdata")}: Class initialization.}
    \item{plot}{\code{signature(x = "HWEdata", y = "missing")}: Provides a pictorial representation for a sample of genotype counts.}
    \item{summary}{\code{signature(object = "HWEdata")}: Extracts the contents of an \code{\link{HWEdata}} object.}
	 }
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\author{Sergio Venturini \email{sergio.venturini@unicatt.it} }
\seealso{
   \code{\link{hwe.ibf}}, 
   \code{\link{hwe.ibf.mc}},
   \code{\link{hwe.ibf.plot}}.
}
\examples{
data.tmp <- c(3, 9, 8)
dataset <- new("HWEdata", data = data.tmp)
}
\keyword{classes}
\keyword{methods}
