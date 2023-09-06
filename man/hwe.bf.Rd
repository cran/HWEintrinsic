\name{hwe.bf}
\alias{hwe.bf}
\title{Standard Bayes Factor for the Hardy-Weinebrg Testing Problem}
\description{
  This function provides the calculations for obtaining the standard Bayes factor for the Hardy-Weinberg testing problem. It implements a common default prior (constant) for both the null and the alternative models.  
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
\note{
  The Bayes factor computed here is for the unrestricted model (\eqn{M_1}{M_1}) against the Hardy-Weinberg case (\eqn{M_0}{M_0}).
}
\references{
  Consonni, G., Moreno, E., and Venturini, S. (2011). "Testing Hardy-Weinberg equilibrium: an objective Bayesian analysis". Statistics in Medicine, \bold{30}, 62--74. \url{https://onlinelibrary.wiley.com/doi/10.1002/sim.4084/abstract}
}
\author{Sergio Venturini \email{sergio.venturini@unicatt.it} }
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
