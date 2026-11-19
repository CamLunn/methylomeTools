\name{mask_matrices_by_coverage}
\alias{mask_matrices_by_coverage}
\title{Methylation utility: mask matrices by coverage.}
\description{Methylation utility: mask matrices by coverage.}
\usage{
mask_matrices_by_coverage(cov_mat,
    beta_mat = NULL,
    meth_mat = NULL,
    extra_mats = NULL,
    lo.count = NULL,
    hi.count = NULL,
    lo.perc = NULL,
    hi.perc = NULL,
    mode = c("mask", "subset"),
    keep_rule = c("any", "all", "max_missing"),
    max_missing_prop = 1.0,
    return_what = c("beta", "meth", "cov", "all"),
    verbose = TRUE)
}
\arguments{
  \item{cov_mat}{Total coverage matrix, with CpGs in rows and samples in columns.}
  \item{beta_mat}{Configuration argument; see the default in the usage and the module guide.}
  \item{meth_mat}{Methylation matrix in SVA functions; methylated read counts in coverage utilities. See module guide.}
  \item{extra_mats}{Configuration argument; see the default in the usage and the module guide.}
  \item{lo.count}{Configuration argument; see the default in the usage and the module guide.}
  \item{hi.count}{Configuration argument; see the default in the usage and the module guide.}
  \item{lo.perc}{Configuration argument; see the default in the usage and the module guide.}
  \item{hi.perc}{Configuration argument; see the default in the usage and the module guide.}
  \item{mode}{Configuration argument; see the default in the usage and the module guide.}
  \item{keep_rule}{Configuration argument; see the default in the usage and the module guide.}
  \item{max_missing_prop}{Configuration argument; see the default in the usage and the module guide.}
  \item{return_what}{Configuration argument; see the default in the usage and the module guide.}
  \item{verbose}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A matrix, list of matrices, vector or summary table depending on the utility.}
\details{See the README and docs directory for input formats and optional dependencies.}
