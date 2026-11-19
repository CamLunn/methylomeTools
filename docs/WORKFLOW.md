# Workflow and input conventions

## Prepare aligned inputs

Matrices have features in rows and samples in columns. Sample IDs and feature IDs
are unique, and metadata follows matrix column order. Missing covariates can cause
model.matrix to omit rows.

```r
metadata <- metadata[match(colnames(beta_matrix), metadata$sample_id), , drop=FALSE]
stopifnot(identical(metadata$sample_id, colnames(beta_matrix)))
rownames(metadata) <- metadata$sample_id
mod <- model.matrix(~ exposure + age + sex + batch, metadata)
mod0 <- model.matrix(~ age + sex + batch, metadata)
stopifnot(nrow(mod) == ncol(beta_matrix), qr(mod)$rank == ncol(mod),
          qr(mod0)$rank == ncol(mod0))
sv <- build_methyl_svs(beta_matrix, mod, mod0, value_type="beta",
  feature_strategy="hvf_only", n_hvf=200000, override_hvf=TRUE,
  sva_method="irw", num_sv_method="be")
```

Pass value_type="M" when supplying already transformed logit values. Feature
preparation can be cached with prepare_methyl_matrix and reused with
run_sva_on_matrix for different models on the same cohort. Generic probe IDs are
accepted; genomic coordinates in chr.pos or chr.start.end form are optional metadata.

Imputation is not implemented: leave impute=FALSE. With that default, features
containing any remaining NA after selection are removed. A missingness_max of
0.05 therefore does not mean incomplete selected features are retained for fitting.
Preprocessing and SVA each run in one R process. The workers argument is retained
for compatibility; values above one warn. Threaded BLAS may accelerate matrix
operations through blas_threads, depending on the installed library. The legacy
two_step option is rejected because it called svaseq rather than a two-step fit.
Use irw for methylation M values.

Feature ranking uses the supplied scale. Selected beta values are converted to
log2 M values with endpoint clipping at 1e-6 before centering or fitting. Older
snapshots did not perform this conversion, so beta-input results differ and older
beta-input cached preparations are rejected. Cached preparations from the corrected version are reusable
only for the same samples and preprocessing settings.

## Join SVs and run DSS

The returned n_sv records the estimated SV count. The following join uses sample
IDs to preserve metadata order.

```r
idx <- match(metadata$sample_id, sv$sv$sample_id)
stopifnot(!anyNA(idx))
sv_names <- grep("^SV[0-9]+$", names(sv$sv), value=TRUE)
metadata[sv_names] <- sv$sv[idx, sv_names, drop=FALSE]
fit_formula <- reformulate(c("exposure", "age", "sex", "batch", sv_names))
run <- run_dss_parallel(bs=bs, meta=metadata, formula=fit_formula,
  outcome_col="exposure", outcome_type="continuous", sample_id_col="sample_id",
  n_workers=4, parallel_backend="psock", smoothing=FALSE,
  partition_mode="chromosome", outdir="outputs/exposure_dss", result_prefix="exposure")
results <- readRDS("outputs/exposure_dss/exposure_merged_results.rds")$main
```

DSS requires methylated and total read counts; array beta values alone are
insufficient. SVA adjustment does not replace an explicit biological model.
canonical_chr_only=TRUE retains chr1..chr22, chrX and chrY by default. For other
species or chromosome conventions, disable it and set chr_subset explicitly.

Global BH correction is performed over merged P values separately for each named
test. Per-unit adjusted P values are retained for reference. The correction
universe is the analyzed merged sites, not all genomic CpGs. Job summaries record
failed units. Block and chromosome partitioning fit separate subsets and can
produce different estimates.
Smoothing forces chromosome partitioning to avoid block boundaries.

Use fresh output directories. resume=TRUE is only for the exact same input counts
and settings; the configuration guard compares metadata, loci and model settings,
but does not hash the complete count matrices. Changed counts require a fresh run.

## Regions and plotting

```r
regions <- run_combp_on_dss(results, outdir="outputs/combp",
  dist_cutoff=300, bin_size=60, seed=0.01, n_cores=4)$dmrs
# Example significance and minimum-site thresholds
regions <- subset(regions, sidak < 0.05 & nprobe >= 3)
sites <- make_dmr_site_table(bs, regions)
long <- make_dmr_methylation_long(bs, sites, dmr_id=sites$dmr_id[1],
  metadata=metadata, min_coverage=10, flank_radius_bp=500)
browser <- plot_dmr_browser(long$dmr_df, long$methylation_long,
  dmr_id=sites$dmr_id[1], phenotype_col="exposure", curve_method="none")
browser$plot
```

Empty comb-p results or zero site overlaps provide no region for browser
selection. plot_top_dml_betas_scatter and plot_top_dml_betas_boxplot accept a BSseq
object plus DML results and metadata. Their extracted beta values are unadjusted
observations, even when the DML statistics came from an adjusted model. DMR curves
are descriptive. Optional browser tracks require RDS annotation tables in matching
genome assembly, chromosome style and coordinate convention. A plotting window
cannot show neighboring data that were not extracted: coordinate padding and
flank_radius_bp should be chosen together.

The hg38 annotation helper remains explicitly human-specific. Other annotation assemblies require compatible external resources.
