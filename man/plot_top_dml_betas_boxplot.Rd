\name{plot_top_dml_betas_boxplot}
\alias{plot_top_dml_betas_boxplot}
\title{Plot selected sequencing CpG beta values by phenotype tertiles.}
\description{Plot selected sequencing CpG beta values by phenotype tertiles.}
\usage{
plot_top_dml_betas_boxplot(bs,
    dml_table,
    meta,
    sample_id_col = "sample_id",
    group_col,
    fdr_col = "fdr_global",
    chr_col = "chr",
    pos_col = "pos",
    top_n = 9,
    cpg_ids = NULL,
    min_cov = 1,
    tertile_labels = c("Low", "Middle", "High"),
    pdf_file = NULL,
    ncol = 3,
    point_alpha = 0.75,
    point_size = 1.8)
}
\arguments{
  \item{bs}{A BSseq object containing methylated counts and total coverage.}
  \item{dml_table}{Configuration argument; see the default in the usage and the module guide.}
  \item{meta}{Sample metadata, with one row per sample.}
  \item{sample_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{group_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{fdr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{pos_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{top_n}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpg_ids}{Configuration argument; see the default in the usage and the module guide.}
  \item{min_cov}{Configuration argument; see the default in the usage and the module guide.}
  \item{tertile_labels}{Configuration argument; see the default in the usage and the module guide.}
  \item{pdf_file}{Configuration argument; see the default in the usage and the module guide.}
  \item{ncol}{Configuration argument; see the default in the usage and the module guide.}
  \item{point_alpha}{Configuration argument; see the default in the usage and the module guide.}
  \item{point_size}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A list containing plot and extracted or summarized data.}
\details{See the README and docs directory for input formats and optional dependencies.}
