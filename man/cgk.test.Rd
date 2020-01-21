\name{cgk.test}
\alias{cgk.test}
\title{
Test of independence based on copula and Gaussian kernel
}
\description{
Runs a test for independence among several univariate random variables based on copula and Gaussian kernel.
}
\usage{
cgk.test(X, method = c("sum", "max", "fdr"), alpha = 0.05, B.parameters = 1000,
B = 1000, parameters = NULL, preset = NULL, test.preset = NULL)
}
\arguments{
  \item{X}{
A matrix whose rows represent observations.
}
  \item{method}{
A vector of strings containting any combination of these three strings: \code{"sum"}, \code{"max"} and \code{"fdr"}. Depending upon those combination, \code{cgk.test} will perform those methods. For more details, see reference.
}
  \item{alpha}{
A numeric value between 0 and 1 specifying the level of the test.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
}
  \item{B}{
An integer value specifying the number of Monte-Carlo iterations to be used for determining the cutoff of the test.
}
  \item{parameters}{
If \code{NULL}, then \code{cgk.test.preset} automatically calculates the parameters. Otherwise it calculates preset using the given set of parameters.
}
  \item{preset}{
If the list of numeric vectors is precalculated by the \code{cgk.preset} function, supply that list here.
}
  \item{test.preset}{
If the list of numeric vectors is precalculated by the \code{cgk.test.preset} function, supply that list here.
}
}
\details{
See reference.
}
\value{
  \item{sum$reject}{
Boolean value indicating whether the null hypothesis should be rejected or not when hypothesis testing is done by the Tsum statistics.
}
  \item{sum$p.value}{
Pvalue of the Tsum statistics.
}
  \item{max$reject}{
Boolean value indicating whether the null hypothesis should be rejected or not when hypothesis testing is done by the Tmax statistics.
}
  \item{max$p.value}{
Pvalue of the Tmax statistics.
}
  \item{fdr$reject}{
Boolean value indicating whether the null hypothesis should be rejected or not when hypothesis testing is done by the Tfdr statistics.
}
}
\references{
Roy et al. (2019)
}
\author{
Angshuman Roy
}
\examples{
set.seed(1)
x <- rnorm(10)
y <- rnorm(10)
z <- x+y+rnorm(10,sd=0.1)
X<- cbind(x,y,z)
set.seed(1)
cgk.test(X,B=100,alpha=0.05)
}
\keyword{nonparametric}
\keyword{htest}
