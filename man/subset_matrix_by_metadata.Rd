\name{subset_matrix_by_metadata}
\alias{subset_matrix_by_metadata}
\title{Methylation utility: subset matrix by metadata.}
\description{Methylation utility: subset matrix by metadata.}
\usage{
subset_matrix_by_metadata(mat,
    strict = TRUE,
    verbose = TRUE,
    maf_pct_max = NULL,
    maf_pct_min = NULL,
    has_snp = NULL,
    cpg_context_include = NULL,
    cpg_context_exclude = NULL,
    epigenetic_include = NULL,
    epigenetic_exclude = NULL,
    genic_include = NULL,
    genic_exclude = NULL,
    gene_symbols_include = NULL,
    gene_symbols_exclude = NULL,
    tss_abs_max = NULL,
    tss_abs_min = NULL,
    return = c("matrix", "metadata", "both"),
    db_file = NULL,
    table_name = NULL)
}
\arguments{
  \item{mat}{Configuration argument; see the default in the usage and the module guide.}
  \item{strict}{Configuration argument; see the default in the usage and the module guide.}
  \item{verbose}{Configuration argument; see the default in the usage and the module guide.}
  \item{maf_pct_max}{Configuration argument; see the default in the usage and the module guide.}
  \item{maf_pct_min}{Configuration argument; see the default in the usage and the module guide.}
  \item{has_snp}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpg_context_include}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpg_context_exclude}{Configuration argument; see the default in the usage and the module guide.}
  \item{epigenetic_include}{Configuration argument; see the default in the usage and the module guide.}
  \item{epigenetic_exclude}{Configuration argument; see the default in the usage and the module guide.}
  \item{genic_include}{Configuration argument; see the default in the usage and the module guide.}
  \item{genic_exclude}{Configuration argument; see the default in the usage and the module guide.}
  \item{gene_symbols_include}{Configuration argument; see the default in the usage and the module guide.}
  \item{gene_symbols_exclude}{Configuration argument; see the default in the usage and the module guide.}
  \item{tss_abs_max}{Configuration argument; see the default in the usage and the module guide.}
  \item{tss_abs_min}{Configuration argument; see the default in the usage and the module guide.}
  \item{return}{Configuration argument; see the default in the usage and the module guide.}
  \item{db_file}{User-provided SQLite annotation database path; no database is bundled.}
  \item{table_name}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A read-only connection, metadata table, filtered matrix or database information.}
\details{See the README and docs directory for input formats and optional dependencies.}
