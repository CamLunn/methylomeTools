# Compatibility aliases preserve existing entry points without duplicating code.
run_dss_multifactor_by_chr <- function(...) run_dss_parallel(...)
make_truseq_mats_by_mincov_from_base <- function(...) make_coverage_threshold_matrices(...)

# Retain the historical array mapping layout for existing callers.
process_beta_EPIC <- function(beta_EPIC_i, EPIC_qc, truseq_epic_map_df,
                              sample_id_col="Basename", probe_col="arraymap_probe_epic_v1") {
  mapping <- truseq_epic_map_df
  mapping$.mapping_id <- as.character(mapping[[probe_col]])
  align_array_beta_to_cpgs(beta_EPIC_i, mapping, probe_col=probe_col,
    cpg_id_col=".mapping_id", sample_metadata=EPIC_qc,
    array_sample_col=sample_id_col, sample_id_col="sample_id", na.rm=FALSE)
}
