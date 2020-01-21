\name{cgk.test.preset}
\alias{cgk.test.preset}
\title{
A list of numeric vectors pre-calculated for \code{cgk.test} function
}
\description{
When \code{cgk.test} function has to be calculated repeatedly, we can save a lot of computation by pre-calculating some numeric vectors that \code{cgk.test} function uses by this \code{cgk.test.preset} function.
}
\usage{
cgk.test.preset(n, d, parameters = NULL, B.parameters = 1000, B = 1000, preset = NULL)
}
\arguments{
  \item{n}{
Sample size.
}
  \item{d}{
Number of univariate random variables.
}
  \item{parameters}{
If \code{NULL}, then \code{cgk.test.preset} automatically calculates the parameters. Otherwise it calculates preset using the given set of parameters.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
}
  \item{B}{
An integer value specifying the number of Monte-Carlo iterations to be used for determining the cutoff of the test.
}
  \item{preset}{
If the list of numeric vectors is precalculated by the \code{cgk.preset} function, supply that list here.
}
}
\details{
See reference.
}
\value{
A list of numeric vectors.
}
\references{
Roy et al. (2019)
}
\author{
Angshuman Roy
}
\examples{
#automated parameters
set.seed(0)
cgk.test.preset(10,2,B.parameters=1000,B=1000)

#manual parameters
cgk.test.preset(10,2,parameters=c(0.3,0.5,1),B=1000)
}
\keyword{nonparametric}
