\name{hwe.ibf.mc}
\alias{hwe.ibf.mc}
\title{Testing Hardy-Weinberg Equilibrium Using an Intrinsic Prior Approach}
\description{
  This function implements the Monte Carlo estimation of the Bayes factor based on intrinsic priors for the Hardy-Weinberg testing problem as described in Consonni et al. (2010).
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
\references{
  Consonni, G., Moreno, E. and Venturini, S. (2010), "Testing Hardy-Weinberg Equilibrium: an Objective Bayesian Analysis". Technical Report. \url{http://economia.unipv.it/~gconsonni/www/papers/HardyWeinberg.pdf}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{hwe.ibf}}, 
   \code{\link{hwe.ibf.plot}}.
}
\examples{
# Example 1 #

data(GuoThompson9)
plot(GuoThompson9)
n <- sum(GuoThompson9@vec, na.rm = TRUE)
res <- hwe.ibf.mc(GuoThompson9, t = n/2, M = 100000, verbose = TRUE)
summary(res, plot = TRUE)

# Example 2 #
\donttest{
# ATTENTION: the following code may take a long time to run! #

M <- 300000
f <- seq(.1, 1, .05)
n <- sum(GuoThompson9@vec, na.rm = TRUE)
out <- hwe.ibf.plot(y = GuoThompson9, t.vec = round(f*n), M = M)
}
}
\keyword{Bayes Factor}
\keyword{intrinsic prior}
\keyword{Hardy-Weinberg}
\keyword{HWE}
