\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY,ANY-method}
\alias{plot,HWEdata,missing-method}
\alias{plot,HWEintr,missing-method}
\title{Plot of Hardy-Weinberg Data and Analysis}
\description{
   Methods for function \code{plot} in Package `graphics' to be used with "\link{HWEdata}" and "\link{HWEintr}" objects.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "HWEdata", y = "missing")}}{
  Pictorial representation for a "\link{HWEdata}" object.
}

\item{\code{signature(x = "HWEintr", y = "missing")}}{
  Graphical representation of Monte Carlo sums for a "\link{HWEintr}" object.
}
}}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{http://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{hwe.ibf}}, 
   \code{\link{hwe.ibf.mc}},
   \code{\link{hwe.ibf.plot}}.
}
\keyword{methods}
\keyword{intrinsic prior}
