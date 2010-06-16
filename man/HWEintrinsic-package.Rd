\name{HWEintrinsic-package}
\alias{HWEintrinsic-package}
\alias{HWEintrinsic}
\docType{package}
\title{Hardy-Weinberg Equilibrium Testing using an Objective Bayesian Approach}
\description{
  This package implements the intrinsic prior methodology to the problem of testing whether a population is in the Hardy-Weinberg equilibrium. The class of intrinsic priors under the full model is indexed by a tuning quantity, the training sample size. The analysis is robust if the decision of rejecting the null hypothesis does not change as the intrinsic prior varies over this class.  
}
\details{
\tabular{ll}{
Package: \tab HWEintrinsic\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-06-16\cr
License: \tab GPL 2 or greater\cr
LazyLoad: \tab yes\cr
}
  The package is loaded with the usual \code{library(HWEintrinsic)} command. The most important functions are \code{\link{hwe.ibf}}, \code{\link{hwe.ibf.mc}}, and \code{\link{hwe.ibf.plot}}.
  It contains also few data sets which have been extensively used in the literature.   
}
\references{
  Consonni, G., Moreno, E. and Venturini, S. (2010), "Testing Hardy-Weinberg Equilibrium: an Objective Bayesian Analysis". Technical Report. \url{http://economia.unipv.it/~gconsonni/www/papers/HardyWeinberg.pdf}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{hwe.ibf}}, 
   \code{\link{hwe.ibf.mc}},
   \code{\link{hwe.ibf.plot}},
   \code{\link{hwe.bf}}.
}
\keyword{intrinsic prior}
\keyword{Bayes Factor}
\keyword{Hardy-Weinberg}
\keyword{HWE}
