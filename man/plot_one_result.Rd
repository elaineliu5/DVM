% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/examples.R
\name{plot_one_result}
\alias{plot_one_result}
\title{plot_one_result}
\usage{
plot_one_result(result, max_range)
}
\arguments{
\item{result}{a returned value from run_simulations}

\item{max_range}{the maximum range of probabilities to be used for setting the color scale (breaks/legends)}
}
\value{
None
}
\description{
plot_one_result plots the logistic probabilities 
and population sample points with a color scale corresponding to 
a given probability range
}
\examples{
plot_one_result(result1, max_range=range(unlist(result$range)))

}
