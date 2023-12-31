% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{average_piks_list}
\alias{average_piks_list}
\title{Average of priors: pnk, pak and pck from list (memory intensive)}
\usage{
average_piks_list(params, nsnps, covar_vec, nits, thin, covar = FALSE)
}
\arguments{
\item{params}{Vector of parameters: \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}}

\item{nsnps}{number of snps}

\item{covar_vec}{Vector of the covariate}

\item{nits}{Number of iterations run in mcmc}

\item{thin}{thinning}

\item{covar}{logical: was the covariate inflormation  used? default: False}
}
\value{
average pik matrix of priors: pnk, pak and pck
}
\description{
Average of priors: pnk, pak and pck from list (memory intensive)
}
