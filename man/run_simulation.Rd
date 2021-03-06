% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/examples.R
\name{run_simulation}
\alias{run_simulation}
\title{run_simulation}
\usage{
run_simulation(species_name, month, N, region, resolution = 10)
}
\arguments{
\item{species_name}{a vector containing genus and species name}

\item{month}{an interger (1-12) indicating the month of interest}

\item{N}{number of species samples to draw in the reference month}

\item{region}{a vector of region names specifying a region of interest 
from national levels down to county levels}

\item{resolution}{the spatial resolution of environmental layers
(10, 5, or 2.5 (arc min))}
}
\value{
log: a raster stack of logistic probability raster of the given region by month

sample: a list of sample points (spatial points) drawn from the given region by month

range: the range of logistic probabilities in the given region across all months
useful for setting a common color scale when plotting the probability rasters for all months

N: number of samples drawn from each month 
(proportionate to the reference month sample size 
by ratio of mean logistic probabilities)
}
\description{
run_simulation takes in a species name, a month, 
the number of samples to draw for that month, 
a region and optionally the environment layer resolution. 
It performs the simulation for the species' population sampling 
in the region across different months.
}
\examples{
region1 <- c("USA", "Pennsylvania")
species_name1 <- c("Ixodes", "ricinus", "Linnaeus")
species_name2 <- c("mus", "musculus")
result1 <- run_simulation(species_name1, month, N, region1, resolution)
result2 <- run_simulation(species_name2, month, N, region1, resolution)

}
