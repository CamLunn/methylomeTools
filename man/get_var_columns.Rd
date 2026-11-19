\name{get_var_columns}
\alias{get_var_columns}
\title{Methylation utility: get var columns.}
\description{Methylation utility: get var columns.}
\usage{
get_var_columns(include_cpg_id = FALSE,
    db_file = NULL,
    table_name = NULL)
}
\arguments{
  \item{include_cpg_id}{Configuration argument; see the default in the usage and the module guide.}
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
  \item{table_name}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
