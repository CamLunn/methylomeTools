\name{align_array_beta_to_cpgs}
\alias{align_array_beta_to_cpgs}
\title{Align comma-separated array probe IDs to CpG IDs, averaging available probes per CpG.}
\description{Align comma-separated array probe IDs to CpG IDs, averaging available probes per CpG.}
\usage{
align_array_beta_to_cpgs(beta_matrix,
    mapping,
    probe_col = "probe_ids",
    cpg_id_col = "cpg_id",
    sample_metadata = NULL,
    array_sample_col = "array_sample_id",
    sample_id_col = "sample_id",
    na.rm = TRUE)
}
\arguments{
  \item{beta_matrix}{Numeric beta matrix with unique probe rownames and sample columns.}
  \item{mapping}{Data frame mapping comma-separated probe IDs to unique CpG IDs.}
  \item{probe_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpg_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{sample_metadata}{Configuration argument; see the default in the usage and the module guide.}
  \item{array_sample_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{sample_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{na.rm}{Whether to ignore missing values when summarizing.}
}
\value{A matrix, list of matrices, vector or summary table depending on the utility.}
\details{See the README and docs directory for input formats and optional dependencies.}
