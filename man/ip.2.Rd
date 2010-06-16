\name{ip.2}
\alias{ip.2}
\title{Graph of the (Unconditional) Intrinsic Prior for a Two Alleles Dataset}
\description{
  This function produces the 3D plot for the (unconditional) intrinsic prior based on a sample of two alleles data as described in Consonni et al. (2010).
}
\usage{
  ip.2(t, k = 30)
}
\arguments{
  \item{t}{training sample size.}
  \item{k}{number of grid points for the alleles frequencies at which the prior is evaluated.}
}
\value{
  No object is returned by this function.
}
\references{
  Consonni, G., Moreno, E. and Venturini, S. (2010), "Testing Hardy-Weinberg Equilibrium: an Objective Bayesian Analysis". Technical Report. \url{http://economia.unipv.it/~gconsonni/www/papers/HardyWeinberg.pdf}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\note{
  This function provides the plot of the (unconditional) intrinsic prior only for two alleles data.
}
\seealso{
  \code{\link{cip.2}},
  \code{\link{plot-methods}}.
}
\examples{
### The following code reproduces Figure 3 in Consonni et al. (2010) ###
t <- 30
ip.2(t, 40)
}
\keyword{intrinsic prior}
\keyword{unconditional intrinsic prior}
\keyword{HWE}
\keyword{plot}
