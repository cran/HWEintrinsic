\name{HWEintr-class}
\Rdversion{1.1}
\docType{class}
\alias{HWEintr-class}
\alias{HWEintr}
\alias{initialize,HWEintr-method}

\title{Class "HWEintr". Result for the Hardy-Weinberg Testing Problem Using the Intrinsic Prior Approach.}
\description{
  This class encapsulates the results of a Bayesian objective analysis via intrinsic priors for the Hardy-Weinberg Testing Problem as described in Consonni et al. (2010).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HWEintr", bf, npp, draws, data)}, but most often as the result of a call to \code{\link{hwe.ibf}} or to \code{\link{hwe.ibf.mc}}.
}
\section{Slots}{
  \describe{
    \item{\code{bf}:}{Object of class \code{"numeric"}; Bayes factor based on intrinsic priors.}
    \item{\code{npp}:}{Object of class \code{"numeric"}; posterior probability of the null Hardy-Weinberg model.}
    \item{\code{draws}:}{Object of class \code{"numeric"}; individual terms of the Monte Carlo sum using importance sampling.}
    \item{\code{data}:}{Object of class \code{"matrix"}; original data in matrix form.}
  }
}
\section{Methods}{
  \describe{
    %\item{initialize}{\code{signature(.Object = "HWEintr")}: Class initialization.}
    \item{plot}{\code{signature(x = "HWEintr", y = "missing")}: .}
    \item{summary}{\code{signature(object = "HWEintr")}: .}
	 }
}
\references{
  Consonni, G., Moreno, E. and Venturini, S. (2010), "Testing Hardy-Weinberg Equilibrium: an Objective Bayesian Analysis". Technical Report. \url{http://economia.unipv.it/~gconsonni/www/papers/HardyWeinberg.pdf}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{hwe.ibf}}, 
   \code{\link{hwe.ibf.mc}},
   \code{\link{hwe.ibf.plot}}.
}
%\examples{
%showClass("HWEintr")
%}
\keyword{classes}
\keyword{methods}
