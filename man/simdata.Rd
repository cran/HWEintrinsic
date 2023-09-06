\name{simdata}
\alias{simdata}
\alias{dataset1}
\alias{dataset2}
\alias{dataset3}
\alias{dataset4}
\docType{data}
\title{Two Alleles Simulated Data}
\description{
  Four different samples of \eqn{n = 20}{n = 20} genotype counts simulated under the Hardy-Weinberg equilibrium model. 
}
\usage{data(simdata)}
\format{
  Four objects of class \code{\link{HWEdata}}.
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\examples{
data(simdata)
summary(dataset1)
plot(dataset1)
summary(dataset2)
plot(dataset2)
summary(dataset3)
plot(dataset3)
summary(dataset4)
plot(dataset4)

# The following code reproduces Figure 4 in Consonni et al. (2011) #
\dontrun{
# ATTENTION: it may take a long time to run! #

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
\keyword{datasets}
\keyword{HWE}
