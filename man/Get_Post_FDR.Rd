% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{Get_Post_FDR}
\alias{Get_Post_FDR}
\title{Calculate the FDR for all genes}
\usage{
Get_Post_FDR(Gene, q_0)
}
\arguments{
\item{Gene}{The vector of gene symbols.}

\item{q_0}{The vector of posterior probabilities of P(S_i=-1) estimated from Gibbs.}
}
\value{
The data frame including FDR (\code{FDR}) for all genes(\code{Gene}) .
}
\description{
Calculate the FDR for all genes
}
