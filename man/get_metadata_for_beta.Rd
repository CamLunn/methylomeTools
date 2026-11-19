\name{get_metadata_for_beta}
\alias{get_metadata_for_beta}
\title{Methylation utility: get metadata for beta.}
\description{Methylation utility: get metadata for beta.}
\usage{
get_metadata_for_beta(mat,
    columns = NULL,
    strict = TRUE,
    db_file = NULL,
    table_name = NULL)
}
\arguments{
  \item{mat}{Configuration argument; see the default in the usage and the module guide.}
  \item{columns}{Configuration argument; see the default in the usage and the module guide.}
  \item{strict}{Configuration argument; see the default in the usage and the module guide.}
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
  \item{table_name}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
