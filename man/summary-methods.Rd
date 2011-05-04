\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{summary,HWEdata-method}
\alias{summary,HWEintr-method}
\title{Summary of Hardy-Weinberg Data and Analysis}
\description{
  Methods for function \code{summary} in Package `base' to be used with "\link{HWEdata}" and "\link{HWEintr}" objects.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "HWEdata")}}{
  Extracts the slots of a "\link{HWEdata}" object.
}

\item{\code{signature(object = "HWEintr")}}{
  Extracts the slots of a "\link{HWEintr}" object.
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
