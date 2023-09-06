\name{hwe.ibf}
\alias{hwe.ibf}
\title{Testing Hardy-Weinberg Equilibrium Using an Intrinsic Prior Approach}
\description{
  This function implements the exact calculation of the Bayes factor based on intrinsic priors for the Hardy-Weinberg testing problem as described in Consonni et al. (2011).
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
  The Bayes factor computed here is for the unrestricted model (\eqn{M_1}{M_1}) against the Hardy-Weinberg case (\eqn{M_0}{M_0}). This function provides the output only for the two alleles case.
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\author{Sergio Venturini \email{sergio.venturini@unicatt.it} }
\seealso{
   \code{\link{hwe.ibf.plot}},
   \code{\link{hwe.ibf.mc}}.
}
\examples{
\dontrun{
# ATTENTION: the following code may take a long time to run! #

data(Lindley)
hwe.ibf.exact <- Vectorize(hwe.ibf, "t")
f <- seq(.05, 1, .05)
n <- sum(dataL1@data.vec, na.rm = TRUE)

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
}
\keyword{Bayes Factor}
\keyword{intrinsic prior}
\keyword{Hardy-Weinberg}
\keyword{HWE}
