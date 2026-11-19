\name{make_coverage_threshold_matrices}
\alias{make_coverage_threshold_matrices}
\title{Generate a series of coverage-filtered methylation matrices.}
\description{Generate a series of coverage-filtered methylation matrices.}
\usage{
make_coverage_threshold_matrices(beta_base = NULL,
    meth_mat = NULL,
    cov_mat,
    thresholds = c(5, 10, 20, 30, 50, 100),
    cut_threshold = c("default", "sections"),
    mode = c("subset", "mask"),
    keep_rule = c("any", "all", "max_missing"),
    max_missing_prop = NULL,
    lo.count = NULL,
    hi.count = NULL,
    lo.perc = NULL,
    hi.perc = NULL,
    return_what = c("beta", "meth", "cov", "all"),
    extra_mats = NULL,
    verbose = TRUE)
}
\arguments{
  \item{beta_base}{Configuration argument; see the default in the usage and the module guide.}
  \item{meth_mat}{Methylation matrix in SVA functions; methylated read counts in coverage utilities. See module guide.}
  \item{cov_mat}{Total coverage matrix, with CpGs in rows and samples in columns.}
  \item{thresholds}{Configuration argument; see the default in the usage and the module guide.}
  \item{cut_threshold}{Configuration argument; see the default in the usage and the module guide.}
  \item{mode}{Configuration argument; see the default in the usage and the module guide.}
  \item{keep_rule}{Configuration argument; see the default in the usage and the module guide.}
  \item{max_missing_prop}{Configuration argument; see the default in the usage and the module guide.}
  \item{lo.count}{Configuration argument; see the default in the usage and the module guide.}
  \item{hi.count}{Configuration argument; see the default in the usage and the module guide.}
  \item{lo.perc}{Configuration argument; see the default in the usage and the module guide.}
  \item{hi.perc}{Configuration argument; see the default in the usage and the module guide.}
  \item{return_what}{Configuration argument; see the default in the usage and the module guide.}
  \item{extra_mats}{Configuration argument; see the default in the usage and the module guide.}
  \item{verbose}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A matrix, list of matrices, vector or summary table depending on the utility.}
\details{See the README and docs directory for input formats and optional dependencies.}
