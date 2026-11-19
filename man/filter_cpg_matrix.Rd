\name{filter_cpg_matrix}
\alias{filter_cpg_matrix}
\title{Methylation utility: filter cpg matrix.}
\description{Methylation utility: filter cpg matrix.}
\usage{
filter_cpg_matrix(mat,
    rule = "ct_ga_no_xy",
    strict = TRUE,
    verbose = TRUE,
    db_file = NULL,
    table_name = NULL)
}
\arguments{
  \item{mat}{Configuration argument; see the default in the usage and the module guide.}
  \item{rule}{Configuration argument; see the default in the usage and the module guide.}
  \item{strict}{Configuration argument; see the default in the usage and the module guide.}
  \item{verbose}{Configuration argument; see the default in the usage and the module guide.}
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
  \item{table_name}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
