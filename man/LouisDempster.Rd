\name{LouisDempster}
\alias{LouisDempster}
\docType{data}
\title{Louis and Dempster (1987) Four Alleles Data}
\description{
  Sample of \eqn{n = 45}{n = 45} genotype counts previously discussed in Louis and Dempster (1987) and Guo and Thompson (1992, Example 1). These data are described in Thomson et al. (1986) and concern the antigen class of 45 French type 1 diabetes patients, with the classes being DR1, DR3, DR4, and Y, a fourth class corresponding to all other antigens.
}
\usage{data(LouisDempster)}
\format{
  An object of class \code{\link{HWEdata}}.
}
\source{
  Louis, E. and Dempster, E. (1987), "An Exact Test for Hardy-Weinberg and Multiple Alleles". Biometrics \bold{Vol. 43}, No. 4, 805--811.
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
  Guo, S.W. and Thompson, E.A. (1992), "Performing the Exact Test of Hardy-Weinberg Proportion for Multiple Alleles". Biometrics, \bold{Vol. 49}, No. 2, 361--372.
  Thomson, G., Klitz, W., Louis, E., Lo, S., Bertrams, L., Baur, M., and Neugebauer, M. (1986), "HLA and IDDM predisposition: New aspects". Genetic Epidemiology, \bold{Vol. 1}, No. 2, 363--368.
}
\examples{
# Example 1 #
\dontrun{
# ATTENTION: the following code may take a long time to run! #


data(LouisDempster)
plot(LouisDempster)
n <- sum(LouisDempster@data.vec, na.rm = TRUE)
out <- hwe.ibf.mc(LouisDempster, t = n/2, M = 100000, verbose = TRUE)
summary(out, plot = TRUE)
}

# Example 2 #
\dontrun{
# ATTENTION: the following code may take a long time to run! #

M <- 300000
f <- seq(.1, 1, .05)
n <- sum(LouisDempster@data.vec, na.rm = TRUE)
out <- hwe.ibf.plot(y = LouisDempster, t.vec = round(f*n), M = M)
}
}
\keyword{datasets}
\keyword{HWE}
