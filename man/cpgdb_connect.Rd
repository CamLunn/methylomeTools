\name{cpgdb_connect}
\alias{cpgdb_connect}
\title{Methylation utility: cpgdb connect.}
\description{Methylation utility: cpgdb connect.}
\usage{
cpgdb_connect(db_file = NULL)
}
\arguments{
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
