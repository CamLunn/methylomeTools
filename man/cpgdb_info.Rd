\name{cpgdb_info}
\alias{cpgdb_info}
\title{Methylation utility: cpgdb info.}
\description{Methylation utility: cpgdb info.}
\usage{
cpgdb_info(db_file = NULL)
}
\arguments{
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
