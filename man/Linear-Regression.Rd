% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Leastsquares.R
\name{Linear Regression}
\alias{Linear Regression}
\alias{linreg}
\title{Linear Regression}
\usage{
linreg (formula,data)
}
\arguments{
\item{formula}{linear expression}

\item{data}{data source}
}
\value{
return the coefficients and coefficient names
}
\description{
calculate coefficients of linear regression
}
\examples{
l1 <- linreg(Petal.Length~Species, data = iris)
linreg(mpg~cyl, data = mtcars)
}
