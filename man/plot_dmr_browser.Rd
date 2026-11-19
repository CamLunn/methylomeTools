\name{plot_dmr_browser}
\alias{plot_dmr_browser}
\title{Plot methylation and optional genomic tracks from long-format data.}
\description{Plot methylation and optional genomic tracks from long-format data.}
\usage{
plot_dmr_browser(dmr_df,
    methylation_long,
    dmr_id,
    phenotype_col,
    regulatory_df = NULL,
    cpg_islands_df = NULL,
    gene_features_df = NULL,
    dmr_id_col = "dmr_id",
    sample_col = "sample_id",
    chr_col = "chr",
    pos_col = "pos",
    beta_col = "beta",
    coverage_col = "coverage",
    regulatory_chr_col = "chr",
    regulatory_start_col = "start",
    regulatory_end_col = "end",
    regulatory_type_col = "type",
    cpg_chr_col = "chr",
    cpg_start_col = "start",
    cpg_end_col = "end",
    cpg_name_col = NULL,
    gene_chr_col = "chr",
    gene_start_col = "start",
    gene_end_col = "end",
    gene_name_col = "gene_name",
    gene_id_col = "gene_id",
    transcript_id_col = "transcript_id",
    gene_feature_col = "feature",
    gene_strand_col = "strand",
    canonical_transcript_id_col = "canonical_transcript_id",
    region = NULL,
    padding = 500,
    phenotype_type = c(
      "auto",
      "categorical",
      "continuous"
    ),
    continuous_display = c(
      "continuous",
      "tertiles"
    ),
    group_colors = NULL,
    continuous_palette = c(
      "#313695",
      "#74ADD1",
      "#FFFFBF",
      "#F46D43",
      "#A50026"
    ),
    plot_means = FALSE,
    connect_means = TRUE,
    mean_linewidth = 0.75,
    curve_method = c(
      "loess",
      "gam",
      "none"
    ),
    curve_span = 0.70,
    curve_alpha = 0.55,
    curve_linewidth = 0.55,
    point_alpha = 0.80,
    point_size_range = c(.5, 3),
    coverage_limits = NULL,
    regulatory_colors = c(
      promoter = "#D73027",
      enhancer = "#FFD92F",
      ctcf = "#80CBE5"
    ),
    cpg_island_color = "#31A354",
    exon_color = "#3B3B3B",
    intron_color = "#696969",
    tss_color = "#7B3294",
    dmr_highlight_color = "#FF8C00",
    dmr_highlight_alpha = 0.13,
    tss_arrow_length_prop = 0.025,
    tss_vertical_nudge = 0.18,
    intron_arrow_spacing = NULL,
    intron_arrow_length_prop = 0.35,
    annotation_arrow_size = 0.12,
    title = NULL,
    show_sample_legend = TRUE,
    max_gene_transcripts = 20,
    relative_track_heights = c(
      methylation = 5.5,
      genes = 2.2,
      cpg = 0.65,
      regulatory = 1
    ),
    base_size = 11)
}
\arguments{
  \item{dmr_df}{Region table. Browser plotting expects chr, start and end with 1-based closed coordinates.}
  \item{methylation_long}{Configuration argument; see the default in the usage and the module guide.}
  \item{dmr_id}{Identifier of the single region to extract or plot.}
  \item{phenotype_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{regulatory_df}{Optional path to an RDS table of regulatory annotations.}
  \item{cpg_islands_df}{Optional path to an RDS table of CpG island annotations.}
  \item{gene_features_df}{Optional path to an RDS table of gene and transcript features.}
  \item{dmr_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{sample_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{pos_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{beta_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{coverage_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{regulatory_chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{regulatory_start_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{regulatory_end_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{regulatory_type_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpg_chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpg_start_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpg_end_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpg_name_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{gene_chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{gene_start_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{gene_end_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{gene_name_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{gene_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{transcript_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{gene_feature_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{gene_strand_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{canonical_transcript_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{region}{Configuration argument; see the default in the usage and the module guide.}
  \item{padding}{Configuration argument; see the default in the usage and the module guide.}
  \item{phenotype_type}{Configuration argument; see the default in the usage and the module guide.}
  \item{continuous_display}{Configuration argument; see the default in the usage and the module guide.}
  \item{group_colors}{Configuration argument; see the default in the usage and the module guide.}
  \item{continuous_palette}{Configuration argument; see the default in the usage and the module guide.}
  \item{plot_means}{Configuration argument; see the default in the usage and the module guide.}
  \item{connect_means}{Configuration argument; see the default in the usage and the module guide.}
  \item{mean_linewidth}{Configuration argument; see the default in the usage and the module guide.}
  \item{curve_method}{Configuration argument; see the default in the usage and the module guide.}
  \item{curve_span}{Configuration argument; see the default in the usage and the module guide.}
  \item{curve_alpha}{Configuration argument; see the default in the usage and the module guide.}
  \item{curve_linewidth}{Configuration argument; see the default in the usage and the module guide.}
  \item{point_alpha}{Configuration argument; see the default in the usage and the module guide.}
  \item{point_size_range}{Configuration argument; see the default in the usage and the module guide.}
  \item{coverage_limits}{Configuration argument; see the default in the usage and the module guide.}
  \item{regulatory_colors}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpg_island_color}{Configuration argument; see the default in the usage and the module guide.}
  \item{exon_color}{Configuration argument; see the default in the usage and the module guide.}
  \item{intron_color}{Configuration argument; see the default in the usage and the module guide.}
  \item{tss_color}{Configuration argument; see the default in the usage and the module guide.}
  \item{dmr_highlight_color}{Configuration argument; see the default in the usage and the module guide.}
  \item{dmr_highlight_alpha}{Configuration argument; see the default in the usage and the module guide.}
  \item{tss_arrow_length_prop}{Configuration argument; see the default in the usage and the module guide.}
  \item{tss_vertical_nudge}{Configuration argument; see the default in the usage and the module guide.}
  \item{intron_arrow_spacing}{Configuration argument; see the default in the usage and the module guide.}
  \item{intron_arrow_length_prop}{Configuration argument; see the default in the usage and the module guide.}
  \item{annotation_arrow_size}{Configuration argument; see the default in the usage and the module guide.}
  \item{title}{Configuration argument; see the default in the usage and the module guide.}
  \item{show_sample_legend}{Configuration argument; see the default in the usage and the module guide.}
  \item{max_gene_transcripts}{Configuration argument; see the default in the usage and the module guide.}
  \item{relative_track_heights}{Configuration argument; see the default in the usage and the module guide.}
  \item{base_size}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A list containing plot and plotted data.}
\details{See the README and docs directory for input formats and optional dependencies.}
