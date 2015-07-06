\name{GuoThompson9}
\alias{GuoThompson9}
\docType{data}
\title{Guo and Thompson (1992) Nine Alleles Data}
\description{
  These data are extracted (Guo and Thompson, 1992) from the Rhesus data in Cavalli-Sforza and Bodmer (1971). They consists of information on 8297 individuals with \eqn{r = 9}{r = 9} groups.
}
\usage{data(GuoThompson9)}
\format{
  An object of class \code{\link{HWEdata}}.
}
\source{
  Guo, S.W. and Thompson, E.A. (1992), "Performing the Exact Test of Hardy-Weinberg Proportion for Multiple Alleles". Biometrics, \bold{Vol. 49}, No. 2, 361--372.
}
\references{
  Cavalli-Sforza, L. and Bodmer, W. (1971), "The Genetics of Human Populations". W.H. Freeman and Company, San Francisco.
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{http://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
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
\keyword{datasets}
\keyword{HWE}
