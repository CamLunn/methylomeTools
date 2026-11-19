\name{get_all_cpg_metadata}
\alias{get_all_cpg_metadata}
\title{Methylation utility: get all cpg metadata.}
\description{Methylation utility: get all cpg metadata.}
\usage{
get_all_cpg_metadata(columns = NULL,
    db_file = NULL,
    table_name = NULL)
}
\arguments{
  \item{columns}{Configuration argument; see the default in the usage and the module guide.}
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
  \item{table_name}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
