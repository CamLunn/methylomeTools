# Column names evaluated inside dplyr and ggplot2 data masks.
utils::globalVariables(c(
  ".center", ".chr_clean", ".chr_group_order", ".chr_index",
  ".chr_length", ".chr_max", ".chr_min", ".chr_number",
  ".chr_raw", ".chromosome_group", ".cpg_order", ".cumulative_pos",
  ".direction", ".dml_chr", ".dml_pos", ".dmr_id",
  ".fdr", ".group_value", ".logp", ".logp_raw",
  ".offset", ".p", ".p_plot", ".plot_end",
  ".plot_label", ".plot_start", ".plot_y", ".pos",
  ".sample_id", ".significant", ".stat", "annotation_id",
  "annotation_type", "bs_row", "chr", "chr_order",
  "cpg_annotation", "cpg_context", "cpg_id", "cpg_label",
  "distance_to_gene_bp", "distance_to_tss_bp", "dmr_id", "dmr_index",
  "end", "end_1based", "enhancer_ids", "enhancer_overlap",
  "gene_features", "gene_id", "group_tertile", "group_tertile_num",
  "nearest_gene_chr", "nearest_gene_end", "nearest_gene_entrez", "nearest_gene_start",
  "nearest_gene_strand", "nearest_gene_symbol", "nearest_gene_tss", "overlapping_annotation_types",
  "overlapping_genes", "overlaps_nearest_gene_body", "p", "pos",
  "position_relative_to_gene", "promoter_overlap", "regulatory_annotations", "signed_distance_to_tss_bp",
  "site_index", "start", "start_1based", "transcript_id",
  "width_bp"
))
