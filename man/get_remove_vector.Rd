\name{get_remove_vector}
\alias{get_remove_vector}
\title{Methylation utility: get remove vector.}
\description{Methylation utility: get remove vector.}
\usage{
get_remove_vector(cpg_ids,
    rule = "ct_ga_no_xy",
    strict = TRUE,
    db_file = NULL,
    table_name = NULL)
}
\arguments{
  \item{cpg_ids}{Configuration argument; see the default in the usage and the module guide.}
  \item{rule}{Configuration argument; see the default in the usage and the module guide.}
  \item{strict}{Configuration argument; see the default in the usage and the module guide.}
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
  \item{table_name}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
