% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_beta.R
\name{MCMC_beta}
\alias{MCMC_beta}
\title{Apply the MCMC to the Beta distribution with a multiplicative log-normal proposal density.}
\usage{
MCMC_beta(x0, N, alph = 4, bet = 2, sigma = 0.125)
}
\arguments{
\item{x0}{Initial guess}

\item{N}{how many iterations for which the MCMC should be ran}
}
\description{
Apply the MCMC to the Beta distribution with a multiplicative log-normal proposal density.
}
