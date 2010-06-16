\name{hwe.ibf}
\alias{hwe.ibf}
\title{Testing Hardy-Weinberg Equilibrium Using an Intrinsic Prior Approach}
\description{
  This function implements the exact calculation of the Bayes factor based on intrinsic priors for the Hardy-Weinberg testing problem as described in Consonni et al. (2010).
}
\usage{
hwe.ibf(y, t)
}
\arguments{
  \item{y}{an object of \link{class} "\link{HWEdata}".}
  \item{t}{training sample size.}
}
\details{
  This function implements the exact formula for the Bayes factor based on intrinsic priors.
}
\value{
  \code{hwe.ibf} returns the value of the Bayes factor based on intrinsic priors.
}
\note{
  This function provides the output only for two alleles case.
}
\references{
  Consonni, G., Moreno, E. and Venturini, S. (2010), "Testing Hardy-Weinberg Equilibrium: an Objective Bayesian Analysis". Technical Report. \url{http://economia.unipv.it/~gconsonni/www/papers/HardyWeinberg.pdf}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{hwe.ibf.plot}},
   \code{\link{hwe.ibf.mc}}.
}
\examples{
data(Lindley)
hwe.ibf.exact <- Vectorize(hwe.ibf, "t")
f <- seq(.05, 1, .05)
n <- sum(dataL1@vec, na.rm = TRUE)

# Dataset 1 #
plot(dataL1)
npp.exact <- 1/(1 + hwe.ibf.exact(round(f*n), y = dataL1))
npp.std <- 1/(1 + hwe.bf(dataL1))
plot(f, npp.exact, type="l", lwd = 2, xlab = "f = t/n",
	ylab = "Null posterior probability")
abline(h = npp.std, col = gray(.5), lty = "longdash")

# Dataset 2 #
plot(dataL2)
npp.exact <- 1/(1 + hwe.ibf.exact(round(f*n), y = dataL2))
npp.std <- 1/(1 + hwe.bf(dataL2))
plot(f, npp.exact, type="l", lwd = 2, xlab = "f = t/n",
	ylab = "Null posterior probability")
abline(h = npp.std, col = gray(.5), lty = "longdash")

# Dataset 3 #
plot(dataL3)
npp.exact <- 1/(1 + hwe.ibf.exact(round(f*n), y = dataL3))
npp.std <- 1/(1 + hwe.bf(dataL3))
plot(f, npp.exact, type="l", lwd = 2, xlab = "f = t/n",
	ylab = "Null posterior probability")
abline(h = npp.std, col = gray(.5), lty = "longdash")

# Dataset 4 #
plot(dataL4)
npp.exact <- 1/(1 + hwe.ibf.exact(round(f*n), y = dataL4))
npp.std <- 1/(1 + hwe.bf(dataL4))
plot(f, npp.exact, type="l", lwd = 2, xlab = "f = t/n",
	ylab = "Null posterior probability")
abline(h = npp.std, col = gray(.5), lty = "longdash")
}
\keyword{Bayes Factor}
\keyword{intrinsic prior}
\keyword{Hardy-Weinberg}
\keyword{HWE}
