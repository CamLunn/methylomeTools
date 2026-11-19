# Function reference

Signatures are generated from this development snapshot. See module guides for input conventions.

## annotate_hg38_dmrs

Annotate 1-based closed regions using human hg38 reference annotations.

```r
annotate_hg38_dmrs(dmr_df,
    chr_col = "chr",
    start_col = "start_1based",
    end_col = "end_1based",
    dmr_id_col = NULL,
    include_intergenic_cpg = TRUE,
    extra_annotations = NULL)
```

## run_dss_multifactor_by_chr

Compatibility alias for run_dss_parallel.

```r
run_dss_multifactor_by_chr(...)
```

## make_truseq_mats_by_mincov_from_base

Compatibility alias for make_coverage_threshold_matrices.

```r
make_truseq_mats_by_mincov_from_base(...)
```

## process_beta_EPIC

Compatibility wrapper for the historical array-to-CpG mapping layout.

```r
process_beta_EPIC(beta_EPIC_i,
    EPIC_qc,
    truseq_epic_map_df,
    sample_id_col="Basename",
    probe_col="arraymap_probe_epic_v1")
```

## cpgdb_defaults

Methylation utility: cpgdb defaults.

```r
cpgdb_defaults()
```

## cpgdb_connect

Methylation utility: cpgdb connect.

```r
cpgdb_connect(db_file = NULL)
```

## fetch_cpg_meta

Methylation utility: fetch cpg meta.

```r
fetch_cpg_meta(cpg_ids,
    columns = NULL,
    strict = TRUE,
    db_file = NULL,
    table_name = NULL)
```

## annotate_cpgs

Methylation utility: annotate cpgs.

```r
annotate_cpgs(cpg_ids,
    columns = NULL,
    strict = TRUE,
    db_file = NULL,
    table_name = NULL)
```

## get_metadata_for_beta

Methylation utility: get metadata for beta.

```r
get_metadata_for_beta(mat,
    columns = NULL,
    strict = TRUE,
    db_file = NULL,
    table_name = NULL)
```

## resolve_remove_rule

Methylation utility: resolve remove rule.

```r
resolve_remove_rule(rule = names(.REMOVE_RULES))
```

## get_remove_vector

Methylation utility: get remove vector.

```r
get_remove_vector(cpg_ids,
    rule = "ct_ga_no_xy",
    strict = TRUE,
    db_file = NULL,
    table_name = NULL)
```

## filter_cpg_matrix

Methylation utility: filter cpg matrix.

```r
filter_cpg_matrix(mat,
    rule = "ct_ga_no_xy",
    strict = TRUE,
    verbose = TRUE,
    db_file = NULL,
    table_name = NULL)
```

## subset_matrix_by_metadata

Methylation utility: subset matrix by metadata.

```r
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
```

## get_cpgs_in_region

Methylation utility: get cpgs in region.

```r
get_cpgs_in_region(chrom,
    start = NULL,
    end = NULL,
    columns = NULL,
    db_file = NULL,
    table_name = NULL)
```

## get_all_cpg_metadata

Methylation utility: get all cpg metadata.

```r
get_all_cpg_metadata(columns = NULL,
    db_file = NULL,
    table_name = NULL)
```

## get_cpg_column_names

Methylation utility: get cpg column names.

```r
get_cpg_column_names(prefix = NULL,
    include_cpg_id = FALSE,
    db_file = NULL,
    table_name = NULL)
```

## get_anno_columns

Methylation utility: get anno columns.

```r
get_anno_columns(include_cpg_id = FALSE,
    db_file = NULL,
    table_name = NULL)
```

## get_var_columns

Methylation utility: get var columns.

```r
get_var_columns(include_cpg_id = FALSE,
    db_file = NULL,
    table_name = NULL)
```

## get_arraymap_columns

Methylation utility: get arraymap columns.

```r
get_arraymap_columns(include_cpg_id = FALSE,
    db_file = NULL,
    table_name = NULL)
```

## cpgdb_info

Methylation utility: cpgdb info.

```r
cpgdb_info(db_file = NULL)
```

## plot_top_dml_betas_boxplot

Plot selected sequencing CpG beta values by phenotype tertiles.

```r
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
```

## plot_top_dml_betas_scatter

Plot selected sequencing CpG beta values against a continuous phenotype.

```r
plot_top_dml_betas_scatter(bs,
    dml_table,
    meta,
    sample_id_col = "sample_id",
    group_col,
    color_col = NULL,
    fdr_col = "fdr_global",
    chr_col = "chr",
    pos_col = "pos",
    top_n = 9,
    cpg_ids = NULL,
    min_cov = 1,
    pdf_file = NULL,
    ncol = 3,
    point_alpha = 0.75,
    point_size = 1.8,
    add_lm = TRUE)
```

## run_combp_on_dss

Call regions from site-level P values using ENmix combp, preserving explicit output coordinate columns.

```r
run_combp_on_dss(dss_res,
    outdir = "combp_results",
    p_col = "pvals",
    dist_cutoff = 1000,
    bin_size = 310,
    seed = 0.01,
    n_cores = 10,
    make_plots = FALSE)
```

## make_dmr_methylation_long

Extract methylation, coverage and optional flanking sites for one DMR.

```r
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
```

## plot_dmr_browser

Plot methylation and optional genomic tracks from long-format data.

```r
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
```

## make_dmr_site_table

Associate 1-based closed DMR intervals with loci in a BSseq object for region plotting.

```r
make_dmr_site_table(bs,
    dmrs,
    chr_col="chr",
    start_col="start_1based",
    end_col="end_1based",
    dmr_id_col=NULL)
```

## run_dss_parallel

Fit DSS multifactor models in genomic work units and apply BH correction across merged results for each test.

```r
run_dss_parallel(bs,
    meta,
    formula,
    outcome_col,
    outcome_type = c("auto", "continuous", "binary", "categorical"),
    ref = NULL,
    contrast_structure = NULL,
    multi_group_default = c("term", "one_vs_rest", "none"),
    coef = NULL,
    term = NULL,
    Contrast = NULL,
    sample_id_col = NULL,
    smoothing = FALSE,
    smoothing.span = 500,
    backend = c("live", "slurm"),
    n_workers = .available_workers(),
    parallel_backend = c("psock", "fork"),
    partition_mode = c("block", "chromosome"),
    target_block_size = 150000L,
    pin_blas_threads_n = 1L,
    chr_subset = NULL,
    canonical_chr_only = TRUE,
    outdir = file.path(getwd(), "dss_multifactor_by_chr"),
    result_prefix = "dss_mf",
    save_rds = TRUE,
    save_fit = FALSE,
    save_compression = c("gzip", "none", "xz"),
    save_test_only = TRUE,
    save_log_memory = TRUE,
    auto_merge = TRUE,
    resume = FALSE,
    slurm_partition = NULL,
    slurm_time = "24:00:00",
    slurm_mem = "16G",
    slurm_cpus_per_task = 1,
    slurm_account = NULL,
    slurm_module_load = character(),
    rscript_path = "Rscript",
    slurm_subset_compress = "gzip",
    verbose = TRUE,
    gc_each_chr = TRUE)
```

## plot_dss_manhattan

Plot site-level significance and optional signed statistics along chromosomes.

```r
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
```

## mask_matrices_by_coverage

Methylation utility: mask matrices by coverage.

```r
mask_matrices_by_coverage(cov_mat,
    beta_mat = NULL,
    meth_mat = NULL,
    extra_mats = NULL,
    lo.count = NULL,
    hi.count = NULL,
    lo.perc = NULL,
    hi.perc = NULL,
    mode = c("mask", "subset"),
    keep_rule = c("any", "all", "max_missing"),
    max_missing_prop = 1.0,
    return_what = c("beta", "meth", "cov", "all"),
    verbose = TRUE)
```

## make_coverage_threshold_matrices

Generate a series of coverage-filtered methylation matrices.

```r
make_coverage_threshold_matrices(beta_base = NULL,
    meth_mat = NULL,
    cov_mat,
    thresholds = c(5, 10, 20, 30, 50, 100),
    cut_threshold = c("default", "sections"),
    mode = c("subset", "mask"),
    keep_rule = c("any", "all", "max_missing"),
    max_missing_prop = NULL,
    lo.count = NULL,
    hi.count = NULL,
    lo.perc = NULL,
    hi.perc = NULL,
    return_what = c("beta", "meth", "cov", "all"),
    extra_mats = NULL,
    verbose = TRUE)
```

## align_array_beta_to_cpgs

Align comma-separated array probe IDs to CpG IDs, averaging available probes per CpG.

```r
align_array_beta_to_cpgs(beta_matrix,
    mapping,
    probe_col = "probe_ids",
    cpg_id_col = "cpg_id",
    sample_metadata = NULL,
    array_sample_col = "array_sample_id",
    sample_id_col = "sample_id",
    na.rm = TRUE)
```

## beta_to_m

Methylation utility: beta to m.

```r
beta_to_m(beta_mat,
    eps = 1e-6)
```

## calc_mean_cpg_beta

Methylation utility: calc mean cpg beta.

```r
calc_mean_cpg_beta(beta_mat,
    na.rm = TRUE,
    min_non_missing = 1)
```

## classify_beta_state

Methylation utility: classify beta state.

```r
classify_beta_state(beta_vec,
    hypo_cutoff = 0.30,
    hyper_cutoff = 0.70,
    ordered_factor = TRUE)
```

## summarize_cpg_beta_states

Methylation utility: summarize cpg beta states.

```r
summarize_cpg_beta_states(beta_mat,
    platform_name = "platform",
    hypo_cutoff = 0.30,
    hyper_cutoff = 0.70,
    min_non_missing = 1)
```

## build_methyl_svs

Prepare methylation features and estimate surrogate variables using sva.

```r
build_methyl_svs(meth_mat,
    mod,
    mod0,
    value_type = c("M", "beta"),
    user_cpg_set = NULL,
    feature_strategy = c("user_set_plus_hvf", "hvf_only", "user_set_only"),
    n_hvf = NULL,
    override_hvf = FALSE,
    hvf_stat = c("mad", "var"),
    missingness_max = 0.05,
    impute = FALSE,
    impute_method = c("knn_chr"),
    knn_k = 10,
    num_sv_method = c("be", "leek", "manual"),
    num_sv_manual = NULL,
    sva_method = c("irw", "two_step"),
    center_features = TRUE,
    scale_features = FALSE,
    use_delayed = TRUE,
    workers = NULL,
    blas_threads = NULL,
    seed = 1,
    log_file = NULL,
    verbose = TRUE,
    diagnostics = NULL,
    diagnostics_dir = NULL,
    return_selected_matrix = FALSE,
    return_feature_table = TRUE,
    return_logs = TRUE,
    write_rds = FALSE,
    rds_file = NULL,
    cache_prep_rds = NULL)
```

## prepare_methyl_matrix

Filter and select methylation features independently of the model, for reuse across analyses.

```r
prepare_methyl_matrix(meth_mat,
    value_type = c("M", "beta"),
    user_cpg_set = NULL,
    feature_strategy = c("user_set_plus_hvf", "hvf_only", "user_set_only"),
    n_hvf = NULL,
    override_hvf = FALSE,
    hvf_stat = c("mad", "var"),
    missingness_max = 0.05,
    impute = FALSE,
    impute_method = c("knn_chr"),
    knn_k = 10,
    center_features = TRUE,
    scale_features = FALSE,
    use_delayed = TRUE,
    workers = NULL,
    blas_threads = NULL,
    seed = 1,
    log_file = NULL,
    verbose = TRUE,
    write_rds = FALSE,
    rds_file = NULL,
    .keep_log_open = FALSE)
```

## run_sva_on_matrix

Estimate surrogate variables from a prepared methylation matrix and full and null designs.

```r
run_sva_on_matrix(prep,
    mod,
    mod0,
    num_sv_method = c("be", "leek", "manual"),
    num_sv_manual = NULL,
    sva_method = c("irw", "two_step"),
    blas_threads = NULL,
    seed = 1,
    diagnostics = NULL,
    diagnostics_dir = NULL,
    return_selected_matrix = FALSE,
    return_feature_table = TRUE,
    return_logs = TRUE,
    write_rds = FALSE,
    rds_file = NULL,
    log_file = NULL,
    verbose = TRUE,
    .reuse_log = FALSE)
```

## submit_methyl_svs_slurm

Write a SLURM driver and job script for SVA, with optional submission.

```r
submit_methyl_svs_slurm(meth_mat,
    mod,
    mod0,
    outdir,
    ...,
    sva_lib = NULL,
    staging_dir = file.path(outdir, "slurm_staging"),
    job_name = "methyl_sva",
    cpus_per_task = 16L,
    mem_gb = NULL,
    time = "12:00:00",
    partition = NULL,
    r_module = NULL,
    extra_modules = character(0),
    openblas_so = NULL,
    blas_threads = NULL,
    account = NULL,
    email = NULL,
    sbatch_bin = NULL,
    submit = TRUE,
    dry_run = FALSE)
```
