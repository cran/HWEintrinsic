\name{lik.multin}
\alias{lik.multin}
\title{Utility Function}
\description{
  This function provides the value of the likelihood function of the full (unrestricted) model, as described in Consonni et al. (2011).
}
\usage{
lik.multin(y, p11, p21)
}
\arguments{
  \item{y}{an object of \link{class} "\link{HWEdata}".}
  \item{p11}{gentotype proportion for the alleles pair \eqn{A_{1}A_{1}}{A1A1}.}
  \item{p21}{gentotype proportion for the alleles pair \eqn{A_{2}A_{1}}{A2A1}.}
}
\details{
  This function has been used to generate the likelihood contours that appear in Figure 4 of Consonni et al. (2011) (see the example below). 
}
\value{
  This function returns the numerical value of the likelihood in correspondence of the argument values.
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{http://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\note{
  This function provides the likelihood function value only for the two alleles case.
}
\seealso{
   \code{\link{hwe.ibf}}, 
   \code{\link{hwe.ibf.mc}},
   \code{\link{hwe.ibf.plot}}.
}
\examples{
# The following code reproduces Figure 4 in Consonni et al. (2011) #
\dontrun{
# ATTENTION: it may take a long time to run! #

data(simdata)
n <- sum(dataset1@data.vec, na.rm = TRUE)
f <- c(.1,.5,1)
t <- round(f*n)
p11 <- p21 <- seq(0,1,length.out=100)
ip <- array(NA,c(length(f),length(p11),length(p21)))
for (k in 1:length(f)) {
	ip[k,,] <- outer(X = p11, Y = p21, FUN = Vectorize(ip.tmp), t[k])
	print(paste(k," / ",length(f),sep=""), quote = FALSE)
}

r <- 2
R <- r*(r + 1)/2
l <- 4
tables <- matrix(NA, nrow = R, ncol = l)
tables[, 1] <- dataset1@data.vec
tables[, 2] <- dataset2@data.vec
tables[, 3] <- dataset3@data.vec
tables[, 4] <- dataset4@data.vec
lik <- array(NA, c(l, length(p11), length(p21)))
M <- 300000
par(mfrow = c(4, 4))
for (k in 1:l) {
	y <- new("HWEdata", data = tables[, k])
	lik[k,,] <- lik.multin(y, p11, p21)
	
	nlev <- 10
	for (q in 1:length(f)) {
		contour(p11, p21, ip[q,,], xlab = expression(p[11]),
				ylab = expression(p[21]), nlevels = nlev, col = gray(0),
				main = "", cex.axis = 1.75, cex.lab = 1.75, labcex = 1.4)
		lines(p11^2, 2*p11*(1 - p11), lty = "longdash", col = gray(0), lwd = 2)
		contour(p11, p21, lik[k,,], nlevels = nlev, add = TRUE,
				col = gray(.7), labcex = 1.2)
		abline(a = 1, b = -1, lty = 3, col = gray(.8))
	}
	hwe.ibf.plot(y = y, t.vec = seq(1,n,1), M = M)
}
}
}
\keyword{Bayes Factor}
\keyword{intrinsic prior}
\keyword{Hardy-Weinberg}
\keyword{HWE}
