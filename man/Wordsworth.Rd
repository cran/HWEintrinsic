\name{Wordsworth}
\alias{Wordsworth}
\docType{data}
\title{Wordsworth et al. (1992) Four Alleles Data}
\description{
  Sample of \eqn{n = 230}{n = 230} genotype counts discussed in Lauretto et al. (2009, Example 3). These data come from a rheumatoid arthritis (RA) study performed by Wordsworth et al. (1992), where two hundred and thirty RA patients were genotyped for the HLA-DR locus. The DR4 allele was subdivided into Dw4, Dw14 and other subtypes. DRX represents all non-DR1, non-Dw4, non-Dw14 alleles.
}
\usage{data(Wordsworth)}
\format{
  An object of class \code{\link{HWEdata}}.
}
\source{
  Lauretto, M.S., Nakano, F., Faria, S.R., Pereira, C.A.B. and Stern, J.M. (2009), "A straightforward multiallelic significance test for the Hardy-Weinberg equilibrium law". Genetics and Molecular Biology, \bold{Vol. 32}, No. 3, 619--625.
}
\references{
  Consonni, G., Moreno, E. and Venturini, S. (2010), "Testing Hardy-Weinberg Equilibrium: an Objective Bayesian Analysis". Technical Report. \url{http://economia.unipv.it/~gconsonni/www/papers/HardyWeinberg.pdf}
  Wordsworth, P., Pile, K.D., Buckley, J.D., Lanchbury, J.S.S., Ollier, B., Lathrop, M. and Bell, J.I. (1992), "HLA heterozygosity contributes to susceptibility to rheumatoid arthritis". American Journal of Human Genetics, \bold{Vol. 51}, No. 3, 585--591.
}
\examples{
# Example 1 #

data(Wordsworth)
plot(Wordsworth)
n <- sum(Wordsworth@vec, na.rm = TRUE)
res <- hwe.ibf.mc(Wordsworth, t = n/2, M = 100000, verbose = TRUE)
summary(res, plot = TRUE)

# Example 2 #
\donttest{
# ATTENTION: the following code may take a long time to run! #

data(Wordsworth)
n <- sum(Wordsworth@vec, na.rm = TRUE)
M <- 300000
f <- seq(.1, 1, .05)
out <- hwe.ibf.plot(y = Wordsworth, t.vec = round(f*n), M = M)
}
}
\keyword{datasets}
\keyword{HWE}
