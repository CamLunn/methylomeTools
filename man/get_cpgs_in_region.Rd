\name{get_cpgs_in_region}
\alias{get_cpgs_in_region}
\title{Methylation utility: get cpgs in region.}
\description{Methylation utility: get cpgs in region.}
\usage{
get_cpgs_in_region(chrom,
    start = NULL,
    end = NULL,
    columns = NULL,
    db_file = NULL,
    table_name = NULL)
}
\arguments{
  \item{chrom}{Configuration argument; see the default in the usage and the module guide.}
  \item{start}{Configuration argument; see the default in the usage and the module guide.}
  \item{end}{Configuration argument; see the default in the usage and the module guide.}
  \item{columns}{Configuration argument; see the default in the usage and the module guide.}
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
  \item{table_name}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
