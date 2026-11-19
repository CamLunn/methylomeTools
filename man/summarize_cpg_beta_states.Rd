\name{summarize_cpg_beta_states}
\alias{summarize_cpg_beta_states}
\title{Methylation utility: summarize cpg beta states.}
\description{Methylation utility: summarize cpg beta states.}
\usage{
summarize_cpg_beta_states(beta_mat,
    platform_name = "platform",
    hypo_cutoff = 0.30,
    hyper_cutoff = 0.70,
    min_non_missing = 1)
}
\arguments{
  \item{beta_mat}{Configuration argument; see the default in the usage and the module guide.}
  \item{platform_name}{Configuration argument; see the default in the usage and the module guide.}
  \item{hypo_cutoff}{Configuration argument; see the default in the usage and the module guide.}
  \item{hyper_cutoff}{Configuration argument; see the default in the usage and the module guide.}
  \item{min_non_missing}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A matrix, list of matrices, vector or summary table depending on the utility.}
\details{See the README and docs directory for input formats and optional dependencies.}
