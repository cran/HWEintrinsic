\name{Lindley}
\alias{Lindley}
\docType{data}
\title{Lindley (1988) Two Alleles Data}
\description{
  Four samples of \eqn{n = 100}{n = 100} genotype counts previously discussed in previously analyzed by Lindley (1988). For the first three sets, the classical "exact" test rejects the null hypothesis of Hardy-Weinberg equilibrium with significance level below 3.4\%, whereas for the last data set the Hardy-Weinberg model is not rejected, its p-value being around 20\%.
}
\usage{data(Lindley)}
\format{
  Four objects of class \code{\link{HWEdata}}.
}
\source{
  Consonni, G., Gutierrez-Pena, E. and Veronese, P. (2008), "Compatible priors for Bayesian model comparison with an application to the Hardy-Weinberg equilibrium model". Test, \bold{Vol. 17}, No. 3, 585--605.  
}
\references{
  Consonni, G., Moreno, E. and Venturini, S. (2010), "Testing Hardy-Weinberg Equilibrium: an Objective Bayesian Analysis". Technical Report. \url{http://economia.unipv.it/~gconsonni/www/papers/HardyWeinberg.pdf}
  Guo, S.W. and Thompson, E.A. (1992), "Performing the Exact Test of Hardy-Weinberg Proportion for Multiple Alleles". Biometrics, \bold{Vol. 49}, No. 2, 361--372.
  Lindley D.V. (1988), "Statistical inference concerning Hardy-Weinberg equilibrium". In: Bernardo, J.M., DeGroot, M.H., Lindley, D.V. and Smith, A.F.M. (eds.), "Bayesian statistics 3". Oxford University Press, 307--326.  
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
\keyword{datasets}
\keyword{HWE}
