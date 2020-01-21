\name{cgk}
\alias{cgk}
\title{
Dependency measure based on copula and Gaussian kernel
}
\description{
Dependency measure that measures amount of association among multiple univariate random variables.
}
\usage{
cgk(X, parameters = NULL, preset = NULL, B.parameters = 1000)
}
\arguments{
  \item{X}{
A matrix whose rows represent observations.
}
  \item{preset}{
If the list of numeric vectors is precalculated by the \code{cgk.preset} function, supply that list here.
}
  \item{parameters}{
If \code{NULL}, then \code{cgk} automatically calculates the parameters. Otherwise it computes the dependency measure using the given set of parameters.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
}
}
\details{
See reference.
}
\value{
A numeric vector containing different values of the dependency measure for each different parameter.
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
X <- cbind(x,y,z)
set.seed(1)
cgk(X)
}
\keyword{nonparametric}
