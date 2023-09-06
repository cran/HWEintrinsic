\name{hwe.ibf.mc}
\alias{hwe.ibf.mc}
\title{Testing Hardy-Weinberg Equilibrium Using an Intrinsic Prior Approach}
\description{
  This function implements the Monte Carlo estimation of the Bayes factor based on intrinsic priors for the Hardy-Weinberg testing problem as described in Consonni et al. (2011).
}
\usage{
hwe.ibf.mc(y, t, M = 10000, verbose = TRUE)
}
\arguments{
  \item{y}{an object of \link{class} "\link{HWEdata}".}
  \item{t}{training sample size.}
  \item{M}{number of Monte Carlo iterations.}
  \item{verbose}{logical; if TRUE the function prints the detailed calculation progress.}
}
\details{
  This function implements a Monte Carlo approximation using importance sampling of the Bayes factor based on intrinsic priors.
}
\value{
  \code{hwe.ibf.mc} returns an object of the class "\link{HWEintr}".
}
\note{
  The Bayes factor computed here is for the unrestricted model (\eqn{M_1}{M_1}) against the Hardy-Weinberg case (\eqn{M_0}{M_0}).
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\author{Sergio Venturini \email{sergio.venturini@unicatt.it} }
\seealso{
   \code{\link{hwe.ibf}}, 
   \code{\link{hwe.ibf.plot}}.
}
\examples{
# Example 1 #
\dontrun{
# ATTENTION: the following code may take a long time to run! #


data(GuoThompson9)
plot(GuoThompson9)
n <- sum(GuoThompson9@data.vec, na.rm = TRUE)
out <- hwe.ibf.mc(GuoThompson9, t = n/2, M = 100000, verbose = TRUE)
summary(out, plot = TRUE)
}

# Example 2 #
\dontrun{
# ATTENTION: the following code may take a long time to run! #

M <- 300000
f <- seq(.1, 1, .05)
n <- sum(GuoThompson9@data.vec, na.rm = TRUE)
out <- hwe.ibf.plot(y = GuoThompson9, t.vec = round(f*n), M = M)
}
}
\keyword{Bayes Factor}
\keyword{intrinsic prior}
\keyword{Hardy-Weinberg}
\keyword{HWE}
