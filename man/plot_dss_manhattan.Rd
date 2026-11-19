\name{plot_dss_manhattan}
\alias{plot_dss_manhattan}
\title{Plot site-level significance and optional signed statistics along chromosomes.}
\description{Plot site-level significance and optional signed statistics along chromosomes.}
\usage{
plot_dss_manhattan(dss_results,
    chr_col = "chr",
    pos_col = "pos",
    p_col = "pvals",
    stat_col = "stat",
    fdr_col = "fdr_global",
    directional = FALSE,
    fdr_threshold = 0.05,
    nominal_p_threshold = NULL,
    point_size = 0.60,
    point_alpha = 0.60,
    chromosome_colors = c("#404040", "#BDBDBD"),
    hypo_color = "#2166AC",
    hyper_color = "#B2182B",
    highlight_fdr = TRUE,
    label_top_n = 0,
    label_col = NULL,
    label_fdr_only = TRUE,
    max_logp = NULL,
    symmetric_y_max = NULL,
    chromosome_gap = 1e6,
    show_chr_boundaries = TRUE,
    chr_boundary_color = "#D9D9D9",
    chr_boundary_linewidth = 0.30,
    include_chr = NULL,
    exclude_chr = NULL,
    title = NULL,
    subtitle = NULL)
}
\arguments{
  \item{dss_results}{Site-level results table with chromosome, position and significance columns.}
  \item{chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{pos_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{p_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{stat_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{fdr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{directional}{Configuration argument; see the default in the usage and the module guide.}
  \item{fdr_threshold}{Configuration argument; see the default in the usage and the module guide.}
  \item{nominal_p_threshold}{Configuration argument; see the default in the usage and the module guide.}
  \item{point_size}{Configuration argument; see the default in the usage and the module guide.}
  \item{point_alpha}{Configuration argument; see the default in the usage and the module guide.}
  \item{chromosome_colors}{Configuration argument; see the default in the usage and the module guide.}
  \item{hypo_color}{Configuration argument; see the default in the usage and the module guide.}
  \item{hyper_color}{Configuration argument; see the default in the usage and the module guide.}
  \item{highlight_fdr}{Configuration argument; see the default in the usage and the module guide.}
  \item{label_top_n}{Configuration argument; see the default in the usage and the module guide.}
  \item{label_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{label_fdr_only}{Configuration argument; see the default in the usage and the module guide.}
  \item{max_logp}{Configuration argument; see the default in the usage and the module guide.}
  \item{symmetric_y_max}{Configuration argument; see the default in the usage and the module guide.}
  \item{chromosome_gap}{Configuration argument; see the default in the usage and the module guide.}
  \item{show_chr_boundaries}{Configuration argument; see the default in the usage and the module guide.}
  \item{chr_boundary_color}{Configuration argument; see the default in the usage and the module guide.}
  \item{chr_boundary_linewidth}{Configuration argument; see the default in the usage and the module guide.}
  \item{include_chr}{Configuration argument; see the default in the usage and the module guide.}
  \item{exclude_chr}{Configuration argument; see the default in the usage and the module guide.}
  \item{title}{Configuration argument; see the default in the usage and the module guide.}
  \item{subtitle}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A list containing plot and diagnostics.}
\details{See the README and docs directory for input formats and optional dependencies.}
