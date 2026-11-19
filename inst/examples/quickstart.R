library(methylomeTools)
set.seed(42)
# Synthetic sequencing data: CpGs in rows, samples in columns.
n_sites <- 800L; n_samples <- 30L
meta <- data.frame(sample_id=paste0("sample", seq_len(n_samples)),
                   exposure=rep(c(0, 1), each=n_samples/2),
                   age=runif(n_samples, 20, 40))
rownames(meta) <- meta$sample_id
chr <- rep(c("chr1", "chr2"), each=n_sites/2)
pos <- rep(seq(1000L, by=20L, length.out=n_sites/2), 2L)
coverage <- matrix(30L, n_sites, n_samples)
prob <- plogis(matrix(rnorm(n_sites), n_sites, n_samples) +
                outer(rep(c(1, 0), c(40, n_sites-40)), meta$exposure * 0.7))
meth <- matrix(rbinom(length(prob), size=30L, prob=as.vector(prob)), n_sites, n_samples)
dimnames(meth) <- dimnames(coverage) <- list(paste(chr, pos, pos, sep="."), meta$sample_id)
bs <- bsseq::BSseq(chr=chr, pos=pos, M=meth, Cov=coverage, sampleNames=meta$sample_id)
beta <- meth / coverage
mod <- model.matrix(~ exposure + age, meta)
mod0 <- model.matrix(~ age, meta)
# One SV keeps the synthetic example small.
sv <- build_methyl_svs(beta, mod, mod0, value_type="beta",
  feature_strategy="hvf_only", n_hvf=400, override_hvf=TRUE,
  num_sv_method="manual", num_sv_manual=1, diagnostics=NULL)
idx <- match(meta$sample_id, sv$sv$sample_id)
stopifnot(!anyNA(idx))
meta$SV1 <- sv$sv$SV1[idx]
# Use a fresh output directory for each analysis.
outdir <- tempfile("methylome_demo_")
fit <- run_dss_parallel(bs=bs, meta=meta, formula=~ exposure + age + SV1,
  outcome_col="exposure", outcome_type="continuous", sample_id_col="sample_id",
  n_workers=1, partition_mode="chromosome", outdir=outdir,
  result_prefix="example", save_compression="none")
results <- readRDS(file.path(outdir, "example_merged_results.rds"))$main
print(plot_dss_manhattan(results)$plot)
print(plot_top_dml_betas_scatter(bs, results, meta,
  group_col="age", top_n=4)$plot)
# An illustrative region for plotting, not a statistically called DMR.
regions <- data.frame(chr="chr1", start_1based=1000L, end_1based=1200L)
sites <- make_dmr_site_table(bs, regions)
long <- make_dmr_methylation_long(bs, sites, dmr_id=sites$dmr_id[1],
  metadata=meta, min_coverage=1, flank_radius_bp=100)
if (requireNamespace("patchwork", quietly=TRUE))
  print(plot_dmr_browser(long$dmr_df, long$methylation_long,
    dmr_id=sites$dmr_id[1], phenotype_col="age", curve_method="none")$plot)
