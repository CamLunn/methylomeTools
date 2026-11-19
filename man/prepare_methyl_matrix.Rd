\name{prepare_methyl_matrix}
\alias{prepare_methyl_matrix}
\title{Filter and select methylation features independently of the model, for reuse across analyses.}
\description{Filter and select methylation features independently of the model, for reuse across analyses.}
\usage{
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
}
\arguments{
  \item{meth_mat}{Methylation matrix in SVA functions; methylated read counts in coverage utilities. See module guide.}
  \item{value_type}{Input scale. Selected beta values are converted to log2 M values with endpoint clipping at 1e-6; M values are used directly. Feature selection uses the supplied scale.}
  \item{user_cpg_set}{Configuration argument; see the default in the usage and the module guide.}
  \item{feature_strategy}{Configuration argument; see the default in the usage and the module guide.}
  \item{n_hvf}{Number of highly variable features to select. Ranking uses the supplied input scale.}
  \item{override_hvf}{Configuration argument; see the default in the usage and the module guide.}
  \item{hvf_stat}{Configuration argument; see the default in the usage and the module guide.}
  \item{missingness_max}{Maximum missing fraction for feature eligibility. With impute=FALSE, selected incomplete features are subsequently removed.}
  \item{impute}{Leave FALSE: chromosome-wise KNN imputation is not implemented.}
  \item{impute_method}{Configuration argument; see the default in the usage and the module guide.}
  \item{knn_k}{Configuration argument; see the default in the usage and the module guide.}
  \item{center_features}{Configuration argument; see the default in the usage and the module guide.}
  \item{scale_features}{Configuration argument; see the default in the usage and the module guide.}
  \item{use_delayed}{Configuration argument; see the default in the usage and the module guide.}
  \item{workers}{Compatibility argument. Preprocessing uses one R process; requests above one warn. Use blas_threads for BLAS threading.}
  \item{blas_threads}{Requested positive integer BLAS thread count. Actual threading depends on the BLAS library; this changes session thread settings. For SLURM, cannot exceed allocated CPUs.}
  \item{seed}{Random seed for reproducible feature diagnostics and SVA estimation.}
  \item{log_file}{Configuration argument; see the default in the usage and the module guide.}
  \item{verbose}{Configuration argument; see the default in the usage and the module guide.}
  \item{write_rds}{Configuration argument; see the default in the usage and the module guide.}
  \item{rds_file}{Configuration argument; see the default in the usage and the module guide.}
  \item{.keep_log_open}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A preparation object, surrogate-variable result list, or SLURM job paths, depending on the function.}
\details{See the README and docs directory for input formats and optional dependencies.}
