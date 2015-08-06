\name{constrKriging-package}
\alias{constrKriging-package}
\alias{constrKriging}
\docType{package}
\title{
Constraint kriging models for various constraints types: monotonicity, boundaries and convexity.
}
\description{
Build constrained kriging models, with prediction, simulation and plot features.

}
\details{
\tabular{ll}{
Package: \tab constrKriging\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-19\cr
License: \tab What license is it under?\cr
}
}
\author{
Hassan Maatouk, Yann Richet, Olivier Roustant
Maintainer: Hassan Maatouk <hassan.maatouk@mines-stetienne.fr>

}
\references{
Hassan Maatouk, Xavier Bay. GAUSSIAN PROCESS EMULATORS FOR COMPUTER EXPERIMENTS WITH INEQUALITY CONSTRAINTS. 2014. <hal-01096751> 
}

\keyword{ package, kriging, constraints, simulation, splines }

\examples{
design = c(0.1, 0.5, 0.9)
response = c(10, 5, 9)
model = kmConvex1D(design, response)
plot(object=model, median=TRUE, quantiles=TRUE, minmax=FALSE, col='red')
}
