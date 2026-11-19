\name{make_dmr_methylation_long}
\alias{make_dmr_methylation_long}
\title{Extract methylation, coverage and optional flanking sites for one DMR.}
\description{Extract methylation, coverage and optional flanking sites for one DMR.}
\usage{
make_dmr_methylation_long(bsobj,
    dmr_dms_long,
    dmr_id,
    metadata = NULL,
    dmr_id_col = "dmr_id",
    dms_chr_col = "dms_chr",
    dms_pos_col = "dms_pos",
    dms_row_col = "dms_row",
    dmr_chr_col = "dmr_chr",
    dmr_start_col = "dmr_start",
    dmr_end_col = "dmr_end",
    dmr_width_col = "dmr_width",
    sample_id_col = "sample_id",
    min_coverage = 10L,
    keep_unobserved = TRUE,
    keep_dms_statistics = TRUE,
    flank_radius_bp = 0L)
}
\arguments{
  \item{bsobj}{A BSseq object containing methylated counts and total coverage.}
  \item{dmr_dms_long}{Region-to-site table, for example from make_dmr_site_table().}
  \item{dmr_id}{Identifier of the single region to extract or plot.}
  \item{metadata}{Optional sample metadata with unique sample IDs.}
  \item{dmr_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dms_chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dms_pos_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dms_row_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dmr_chr_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dmr_start_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dmr_end_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{dmr_width_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{sample_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{min_coverage}{Configuration argument; see the default in the usage and the module guide.}
  \item{keep_unobserved}{Configuration argument; see the default in the usage and the module guide.}
  \item{keep_dms_statistics}{Configuration argument; see the default in the usage and the module guide.}
  \item{flank_radius_bp}{Nonnegative bases extracted on each side of the DMR.}
}
\value{A list containing region results or extracted methylation and site tables.}
\details{See the README and docs directory for input formats and optional dependencies.}
