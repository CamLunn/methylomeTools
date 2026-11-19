\name{make_dmr_site_table}
\alias{make_dmr_site_table}
\title{Associate 1-based closed DMR intervals with loci in a BSseq object for region plotting.}
\description{Associate 1-based closed DMR intervals with loci in a BSseq object for region plotting.}
\usage{
make_dmr_site_table(bs,
    dmrs,
    chr_col="chr",
    start_col="start_1based",
    end_col="end_1based",
    dmr_id_col=NULL)
}
\arguments{
  \item{bs}{A BSseq object containing methylated counts and total coverage.}
  \item{dmrs}{Configuration argument; see the default in the usage and the module guide.}
  \item{chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{start_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{end_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dmr_id_col}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A data frame with one row per region-site association.}
\details{See the README and docs directory for input formats and optional dependencies.}
