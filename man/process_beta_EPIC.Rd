\name{process_beta_EPIC}
\alias{process_beta_EPIC}
\title{Compatibility wrapper for the historical array-to-CpG mapping layout.}
\description{Compatibility wrapper for the historical array-to-CpG mapping layout.}
\usage{
process_beta_EPIC(beta_EPIC_i,
    EPIC_qc,
    truseq_epic_map_df,
    sample_id_col="Basename",
    probe_col="arraymap_probe_epic_v1")
}
\arguments{
  \item{beta_EPIC_i}{Configuration argument; see the default in the usage and the module guide.}
  \item{EPIC_qc}{Configuration argument; see the default in the usage and the module guide.}
  \item{truseq_epic_map_df}{Configuration argument; see the default in the usage and the module guide.}
  \item{sample_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{probe_col}{Configuration argument; see the default in the usage and the module guide.}
}
\value{The corresponding implementation result.}
\details{See the README and docs directory for input formats and optional dependencies.}
