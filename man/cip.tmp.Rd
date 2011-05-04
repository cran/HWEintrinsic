\name{cip.tmp}
\alias{cip.tmp}
\title{Utility Function}
\description{
  This function provides the calculations needed to plot the conditional intrinsic prior based for a two alleles dataset as described in Consonni et al. (2011). 
}
\usage{
  cip.tmp(p11, p21, t, p)
}
\arguments{
  \item{p11}{gentotype proportion for the alleles pair \eqn{A_{1}A_{1}}{A1A1}.}
  \item{p21}{gentotype proportion for the alleles pair \eqn{A_{2}A_{1}}{A2A1}.}
  \item{t}{training sample size.}
  \item{p}{allele frequency used to condition the prior upon.}
}
\value{
  \code{cip.tmp} returns the value of the conditional intrinsic prior evaluated at the arguments values.
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{http://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
  \code{\link{cip.2}},
  \code{\link{plot-methods}}.
}
\keyword{intrinsic prior}
\keyword{conditional intrinsic prior}
\keyword{HWE}
\keyword{plot}
