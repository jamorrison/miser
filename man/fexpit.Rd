% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fexpit.R
\name{fexpit}
\alias{fexpit}
\title{Helper function: expanded expit}
\usage{
fexpit(x, sqz = 1e-06)
}
\arguments{
\item{x}{a vector of values between -Inf and +Inf}

\item{sqz}{the amount by which to 'squeeze', default is .000001}
}
\value{
\preformatted{   a vector of values between 0 and 1 inclusive
}
}
\description{
Helper function: expanded expit
}
