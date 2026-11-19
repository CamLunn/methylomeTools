# Build the coordinate adapter required by make_dmr_methylation_long().
make_dmr_site_table <- function(bs, dmrs, chr_col="chr",
                                start_col="start_1based", end_col="end_1based",
                                dmr_id_col=NULL) {
  if (!inherits(bs, "BSseq")) stop("bs must be a BSseq object")
  required <- c(chr_col, start_col, end_col)
  if (!all(required %in% names(dmrs))) stop("Missing DMR coordinate columns")
  chr <- as.character(dmrs[[chr_col]])
  start <- as.integer(dmrs[[start_col]])
  end <- as.integer(dmrs[[end_col]])
  if (anyNA(chr) || anyNA(start) || anyNA(end) || any(start < 1L) || any(start > end))
    stop("DMR coordinates must be valid 1-based closed intervals")
  if (!is.null(dmr_id_col) && !dmr_id_col %in% names(dmrs)) stop("Missing DMR ID column")
  ids <- if (is.null(dmr_id_col)) paste0(chr, ":", start, "-", end) else as.character(dmrs[[dmr_id_col]])
  if (anyNA(ids) || anyDuplicated(ids)) stop("DMR IDs must be unique and nonmissing")
  regions <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
  loci <- SummarizedExperiment::rowRanges(bs)
  hits <- GenomicRanges::findOverlaps(regions, loci, ignore.strand=TRUE)
  i <- S4Vectors::queryHits(hits); j <- S4Vectors::subjectHits(hits)
  data.frame(dmr_id=ids[i], dmr_chr=chr[i], dmr_start=start[i], dmr_end=end[i],
    dmr_width=end[i]-start[i]+1L, dms_chr=as.character(GenomicRanges::seqnames(loci))[j],
    dms_pos=BiocGenerics::start(loci)[j], dms_row=j, stringsAsFactors=FALSE)
}
