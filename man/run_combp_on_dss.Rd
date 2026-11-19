\name{run_combp_on_dss}
\alias{run_combp_on_dss}
\title{Call regions from site-level P values using ENmix combp, preserving explicit output coordinate columns.}
\description{Call regions from site-level P values using ENmix combp, preserving explicit output coordinate columns.}
\usage{
run_combp_on_dss(dss_res,
    outdir = "combp_results",
    p_col = "pvals",
    dist_cutoff = 1000,
    bin_size = 310,
    seed = 0.01,
    n_cores = 10,
    make_plots = FALSE)
}
\arguments{
  \item{dss_res}{Site-level data frame with chr, pos and the selected P value column.}
  \item{outdir}{Configuration argument; see the default in the usage and the module guide.}
  \item{p_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dist_cutoff}{Configuration argument; see the default in the usage and the module guide.}
  \item{bin_size}{Configuration argument; see the default in the usage and the module guide.}
  \item{seed}{Configuration argument; see the default in the usage and the module guide.}
  \item{n_cores}{Configuration argument; see the default in the usage and the module guide.}
  \item{make_plots}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A list containing region results or extracted methylation and site tables.}
\details{See the README and docs directory for input formats and optional dependencies.}
