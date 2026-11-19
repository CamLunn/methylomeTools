\name{fetch_cpg_meta}
\alias{fetch_cpg_meta}
\title{Methylation utility: fetch cpg meta.}
\description{Methylation utility: fetch cpg meta.}
\usage{
fetch_cpg_meta(cpg_ids,
    columns = NULL,
    strict = TRUE,
    db_file = NULL,
    table_name = NULL)
}
\arguments{
  \item{cpg_ids}{Configuration argument; see the default in the usage and the module guide.}
  \item{columns}{Configuration argument; see the default in the usage and the module guide.}
  \item{strict}{Configuration argument; see the default in the usage and the module guide.}
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
  \item{table_name}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
