\name{submit_methyl_svs_slurm}
\alias{submit_methyl_svs_slurm}
\title{Write a SLURM driver and job script for SVA, with optional submission.}
\description{Write a SLURM driver and job script for SVA, with optional submission.}
\usage{
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
}
\arguments{
  \item{meth_mat}{Methylation matrix in SVA functions; methylated read counts in coverage utilities. See module guide.}
  \item{mod}{Full design matrix, with rows in methylation sample order.}
  \item{mod0}{Nested null design matrix, with rows in the same sample order.}
  \item{outdir}{Configuration argument; see the default in the usage and the module guide.}
  \item{...}{Arguments forwarded to the corresponding implementation.}
  \item{sva_lib}{Optional legacy standalone script path. NULL loads the installed package on the compute node.}
  \item{staging_dir}{Configuration argument; see the default in the usage and the module guide.}
  \item{job_name}{Configuration argument; see the default in the usage and the module guide.}
  \item{cpus_per_task}{Configuration argument; see the default in the usage and the module guide.}
  \item{mem_gb}{Configuration argument; see the default in the usage and the module guide.}
  \item{time}{Configuration argument; see the default in the usage and the module guide.}
  \item{partition}{Optional SLURM partition; NULL uses the cluster default.}
  \item{r_module}{Optional cluster-specific R module. NULL uses Rscript from PATH.}
  \item{extra_modules}{Configuration argument; see the default in the usage and the module guide.}
  \item{openblas_so}{Configuration argument; see the default in the usage and the module guide.}
  \item{blas_threads}{Requested positive integer BLAS thread count. Actual threading depends on the BLAS library; this changes session thread settings. For SLURM, cannot exceed allocated CPUs.}
  \item{account}{Configuration argument; see the default in the usage and the module guide.}
  \item{email}{Configuration argument; see the default in the usage and the module guide.}
  \item{sbatch_bin}{Configuration argument; see the default in the usage and the module guide.}
  \item{submit}{Configuration argument; see the default in the usage and the module guide.}
  \item{dry_run}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A preparation object, surrogate-variable result list, or SLURM job paths, depending on the function.}
\details{See the README and docs directory for input formats and optional dependencies.}
