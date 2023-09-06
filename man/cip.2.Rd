\name{cip.2}
\alias{cip.2}
\title{Graph of the Conditional Intrinsic Prior for a Two Alleles Dataset}
\description{
  This function produces the 3D plot for the conditional intrinsic prior based on a sample of two alleles data as described in Consonni et al. (2011).
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
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\author{Sergio Venturini \email{sergio.venturini@unicatt.it} }
\note{
  This function provides the plot of the conditional intrinsic prior only for two alleles data.
}
\seealso{
  \code{\link{ip.2}},
  \code{\link{plot-methods}}.
}
\examples{
\dontrun{
# ATTENTION: the following code may take a long time to run! #

### The following code reproduces Figure 1 in Consonni et al. (2011) ###
p <- 0.5
t <- 5
cip.2(t, p, 30)
t <- 30
cip.2(t, p, 50)
}
}
\keyword{intrinsic prior}
\keyword{conditional intrinsic prior}
\keyword{HWE}
\keyword{plot}
