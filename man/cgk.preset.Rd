\name{cgk.preset}
\alias{cgk.preset}
\title{
A list of numeric vectors pre-calculated for \code{cgk} function
}
\description{
When \code{cgk} function has to be calculated repeatedly, we can save a lot of computation by pre-calculating some numeric vectors that \code{cgk} function uses by this \code{cgk.preset} function.
}
\usage{
cgk.preset(n,d,parameters=NULL,B.parameters=1000)
}
\arguments{
  \item{n}{
Sample size.
}
  \item{d}{
Number of univariate random variables.
}
  \item{parameters}{
If \code{NULL}, then \code{cgk.preset} automatically calculates the parameters. Otherwise it calculates preset using the given set of parameters.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
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
cgk.preset(10,2,B.parameters=1000)

#manual parameters
cgk.preset(10,2,parameters=c(0.3,0.5,1))
}
\keyword{nonparametric}
