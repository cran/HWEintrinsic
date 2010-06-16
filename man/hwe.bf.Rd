\name{hwe.bf}
\alias{hwe.bf}
\title{Standard Bayes Factor for the Hardy-Weinebrg Testing Problem}
\description{
  This function provides the calculations for obtaining the standard Bayes Factor for the Hardy-Weinberg testing problem. It implements a common default prior (constant) for both the null and the alternative models.  
}
\usage{
  hwe.bf(y)
}
\arguments{
  \item{y}{an object of \link{class} "\link{HWEdata}".}
}
\value{
  \code{hwe.bf} returns the standard Bayes Factor value for the Hardy-Weinberg testing problem (see the references for the details).
}
\references{
  Consonni, G., Moreno, E. and Venturini, S. (2010), "Testing Hardy-Weinberg Equilibrium: an Objective Bayesian Analysis". Technical Report. \url{http://economia.unipv.it/~gconsonni/www/papers/HardyWeinberg.pdf}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{hwe.ibf}}, 
   \code{\link{hwe.ibf.mc}},
   \code{\link{hwe.ibf.plot}}.
}
\examples{
# Example 1 #
data(GuoThompson8)
plot(GuoThompson8)
hwe.bf(GuoThompson8)

# Example 2 #
data(LouisDempster)
plot(LouisDempster)
hwe.bf(LouisDempster)
}
\keyword{Bayes Factor}
\keyword{Hardy-Weinberg}
\keyword{HWE}
