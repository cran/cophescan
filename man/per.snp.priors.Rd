% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/singlevar.R
\name{per.snp.priors}
\alias{per.snp.priors}
\title{per.snp.priors}
\usage{
per.snp.priors(
  nsnps,
  p1 = 1e-04,
  p2 = 1e-04,
  p12 = 1e-05,
  pa = NULL,
  pc = NULL
)
}
\arguments{
\item{nsnps}{number of SNPs}

\item{p1}{prior probability a SNP is associated with trait 1, default 1e-4 (coloc prior)}

\item{p2}{prior probability a SNP is associated with trait 2, default 1e-4 (coloc prior)}

\item{p12}{prior probability a SNP is associated with both traits, default 1e-5 (coloc prior)}

\item{pa}{prior probability that a non-query variant is causally associated with the query trait}

\item{pc}{prior probability that the query variant is causally associated with the query trait}
}
\value{
priors at the query variant
}
\description{
Estimate per snp priors
}
\author{
Ichcha Manipur
}
