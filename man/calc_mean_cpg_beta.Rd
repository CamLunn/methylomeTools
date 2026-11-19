\name{calc_mean_cpg_beta}
\alias{calc_mean_cpg_beta}
\title{Methylation utility: calc mean cpg beta.}
\description{Methylation utility: calc mean cpg beta.}
\usage{
calc_mean_cpg_beta(beta_mat,
    na.rm = TRUE,
    min_non_missing = 1)
}
\arguments{
  \item{beta_mat}{Configuration argument; see the default in the usage and the module guide.}
  \item{na.rm}{Whether to ignore missing values when summarizing.}
  \item{min_non_missing}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A matrix, list of matrices, vector or summary table depending on the utility.}
\details{See the README and docs directory for input formats and optional dependencies.}
