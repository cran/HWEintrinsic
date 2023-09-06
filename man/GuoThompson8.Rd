\name{GuoThompson8}
\alias{GuoThompson8}
\docType{data}
\title{Guo and Thompson (1992) Eight Alleles Simulated Data}
\description{
  A sample of size \eqn{n = 30}{n = 30} of genotype frequencies from a population of \eqn{r = 8}{r = 8} alleles simulated under the Hardy-Weinberg equilibrium when the underlying gene frequencies are (.2, .2, .2, .2, .05, .05, .05, .05); see Guo and Thompson (1992), Example 2.
}
\usage{data(GuoThompson8)}
\format{
  An object of class \code{\link{HWEdata}}.
}
\source{
  Guo, S.W. and Thompson, E.A. (1992), "Performing the Exact Test of Hardy-Weinberg Proportion for Multiple Alleles". Biometrics, \bold{Vol. 48}, No. 2, 361--372.
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\examples{
# Example 1 #
\dontrun{
# ATTENTION: the following code may take a long time to run! #


data(GuoThompson8)
plot(GuoThompson8)
n <- sum(GuoThompson8@data.vec, na.rm = TRUE)
out <- hwe.ibf.mc(GuoThompson8, t = n/2, M = 100000, verbose = TRUE)
summary(out, plot = TRUE)
}

# Example 2 #
\dontrun{
# ATTENTION: the following code may take a long time to run! #

M <- 300000
f <- seq(.1, 1, .05)
n <- sum(GuoThompson8@data.vec, na.rm = TRUE)
out <- hwe.ibf.plot(y = GuoThompson8, t.vec = round(f*n), M = M)
}
}
\keyword{datasets}
\keyword{HWE}
