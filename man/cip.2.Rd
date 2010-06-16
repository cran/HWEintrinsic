\name{cip.2}
\alias{cip.2}
\title{Graph of the Conditional Intrinsic Prior for a Two Alleles Dataset}
\description{
  This function produces the 3D plot for the conditional intrinsic prior based on a sample of two alleles data as described in Consonni et al. (2010).
}
\usage{
  cip.2(t, p, k = 30)
}
\arguments{
  \item{t}{training sample size.}
  \item{p}{allele frequency used to condition the prior upon.}
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
  This function provides the plot of the conditional intrinsic prior only for two alleles data.
}
\seealso{
  \code{\link{ip.2}},
  \code{\link{plot-methods}}.
}
\examples{
### The following code reproduces Figure 1 in Consonni et al. (2010) ###
p <- 0.5
t <- 5
cip.2(t, p, 30)
t <- 30
cip.2(t, p, 50)
}
\keyword{intrinsic prior}
\keyword{conditional intrinsic prior}
\keyword{HWE}
\keyword{plot}
