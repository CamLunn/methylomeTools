\name{run_dss_parallel}
\alias{run_dss_parallel}
\title{Fit DSS multifactor models in genomic work units and apply BH correction across merged results for each test.}
\description{Fit DSS multifactor models in genomic work units and apply BH correction across merged results for each test.}
\usage{
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
}
\arguments{
  \item{bs}{A BSseq object containing methylated counts and total coverage.}
  \item{meta}{Sample metadata, with one row per sample.}
  \item{formula}{R model formula, including the tested variable and adjustment terms.}
  \item{outcome_col}{Metadata column identifying the variable to test.}
  \item{outcome_type}{Configuration argument; see the default in the usage and the module guide.}
  \item{ref}{Configuration argument; see the default in the usage and the module guide.}
  \item{contrast_structure}{Configuration argument; see the default in the usage and the module guide.}
  \item{multi_group_default}{Configuration argument; see the default in the usage and the module guide.}
  \item{coef}{Configuration argument; see the default in the usage and the module guide.}
  \item{term}{Configuration argument; see the default in the usage and the module guide.}
  \item{Contrast}{Configuration argument; see the default in the usage and the module guide.}
  \item{sample_id_col}{Configuration argument; see the default in the usage and the module guide.}
  \item{smoothing}{Configuration argument; see the default in the usage and the module guide.}
  \item{smoothing.span}{Configuration argument; see the default in the usage and the module guide.}
  \item{backend}{Configuration argument; see the default in the usage and the module guide.}
  \item{n_workers}{Number of R worker processes, within allocated resources.}
  \item{parallel_backend}{PSOCK by default; fork must be explicitly requested on Unix.}
  \item{partition_mode}{Block or chromosome work units. Smoothing forces chromosome mode.}
  \item{target_block_size}{Configuration argument; see the default in the usage and the module guide.}
  \item{pin_blas_threads_n}{Configuration argument; see the default in the usage and the module guide.}
  \item{chr_subset}{Configuration argument; see the default in the usage and the module guide.}
  \item{canonical_chr_only}{Configuration argument; see the default in the usage and the module guide.}
  \item{outdir}{Configuration argument; see the default in the usage and the module guide.}
  \item{result_prefix}{Configuration argument; see the default in the usage and the module guide.}
  \item{save_rds}{Configuration argument; see the default in the usage and the module guide.}
  \item{save_fit}{Configuration argument; see the default in the usage and the module guide.}
  \item{save_compression}{Configuration argument; see the default in the usage and the module guide.}
  \item{save_test_only}{Configuration argument; see the default in the usage and the module guide.}
  \item{save_log_memory}{Configuration argument; see the default in the usage and the module guide.}
  \item{auto_merge}{Configuration argument; see the default in the usage and the module guide.}
  \item{resume}{Explicitly reuse completed units only when all inputs, including counts, are unchanged.}
  \item{slurm_partition}{Optional SLURM partition; NULL uses the cluster default.}
  \item{slurm_time}{Configuration argument; see the default in the usage and the module guide.}
  \item{slurm_mem}{Configuration argument; see the default in the usage and the module guide.}
  \item{slurm_cpus_per_task}{Configuration argument; see the default in the usage and the module guide.}
  \item{slurm_account}{Configuration argument; see the default in the usage and the module guide.}
  \item{slurm_module_load}{Optional shell setup lines for the cluster. No module is loaded by default.}
  \item{rscript_path}{Configuration argument; see the default in the usage and the module guide.}
  \item{slurm_subset_compress}{Configuration argument; see the default in the usage and the module guide.}
  \item{verbose}{Configuration argument; see the default in the usage and the module guide.}
  \item{gc_each_chr}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A backend summary; merged result tables are written under outdir when merging succeeds.}
\details{See the README and docs directory for input formats and optional dependencies.}
