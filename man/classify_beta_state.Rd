\name{classify_beta_state}
\alias{classify_beta_state}
\title{Methylation utility: classify beta state.}
\description{Methylation utility: classify beta state.}
\usage{
classify_beta_state(beta_vec,
    hypo_cutoff = 0.30,
    hyper_cutoff = 0.70,
    ordered_factor = TRUE)
}
\arguments{
  \item{beta_vec}{Configuration argument; see the default in the usage and the module guide.}
  \item{hypo_cutoff}{Configuration argument; see the default in the usage and the module guide.}
  \item{hyper_cutoff}{Configuration argument; see the default in the usage and the module guide.}
  \item{ordered_factor}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A matrix, list of matrices, vector or summary table depending on the utility.}
\details{See the README and docs directory for input formats and optional dependencies.}
