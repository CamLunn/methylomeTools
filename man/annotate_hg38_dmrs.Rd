\name{annotate_hg38_dmrs}
\alias{annotate_hg38_dmrs}
\title{Annotate 1-based closed regions using human hg38 reference annotations.}
\description{Annotate 1-based closed regions using human hg38 reference annotations.}
\usage{
annotate_hg38_dmrs(dmr_df,
    chr_col = "chr",
    start_col = "start_1based",
    end_col = "end_1based",
    dmr_id_col = NULL,
    include_intergenic_cpg = TRUE,
    extra_annotations = NULL)
}
\arguments{
  \item{dmr_df}{Region table. Browser plotting expects chr, start and end with 1-based closed coordinates.}
  \item{chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{start_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{end_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dmr_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{include_intergenic_cpg}{Configuration argument; see the default in the usage and the module guide.}
  \item{extra_annotations}{Configuration argument; see the default in the usage and the module guide.}
}
\value{The input region table augmented with annotation columns.}
\details{See the README and docs directory for input formats and optional dependencies.}
