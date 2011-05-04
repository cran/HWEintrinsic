\name{ip.tmp}
\alias{ip.tmp}
\title{Utility Function}
\description{
  This function provides the calculations needed to plot the (unconditional) intrinsic prior based for a two alleles dataset as described in Consonni et al. (2011). 
}
\usage{
  ip.tmp(p11, p21, t)
}
\arguments{
  \item{p11}{gentotype proportion for the the pair of alleles \eqn{A_{1}A_{1}}{A1A1}.}
  \item{p21}{gentotype proportion for the the pair of alleles \eqn{A_{2}A_{1}}{A2A1}.}
  \item{t}{training sample size.}
}
\value{
  \code{ip.tmp} returns the value of the (unconditional) intrinsic prior evaluated at the arguments values.
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{http://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
  \code{\link{ip.2}},
  \code{\link{plot-methods}}.
}
\keyword{intrinsic prior}
\keyword{unconditional intrinsic prior}
\keyword{HWE}
\keyword{plot}
