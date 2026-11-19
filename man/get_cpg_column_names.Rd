\name{get_cpg_column_names}
\alias{get_cpg_column_names}
\title{Methylation utility: get cpg column names.}
\description{Methylation utility: get cpg column names.}
\usage{
get_cpg_column_names(prefix = NULL,
    include_cpg_id = FALSE,
    db_file = NULL,
    table_name = NULL)
}
\arguments{
  \item{prefix}{Configuration argument; see the default in the usage and the module guide.}
  \item{include_cpg_id}{Configuration argument; see the default in the usage and the module guide.}
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
  \item{table_name}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
