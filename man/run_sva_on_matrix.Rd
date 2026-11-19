\name{run_sva_on_matrix}
\alias{run_sva_on_matrix}
\title{Estimate surrogate variables from a prepared methylation matrix and full and null designs.}
\description{Estimate surrogate variables from a prepared methylation matrix and full and null designs.}
\usage{
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
}
\arguments{
  \item{prep}{Configuration argument; see the default in the usage and the module guide.}
  \item{mod}{Full design matrix, with rows in methylation sample order.}
  \item{mod0}{Nested null design matrix, with rows in the same sample order.}
  \item{num_sv_method}{Estimate the SV count using be or leek, or supply a manual count.}
  \item{num_sv_manual}{Nonnegative integer SV count when num_sv_method is manual. Must leave residual degrees of freedom.}
  \item{sva_method}{Use irw for sva::sva on methylation M values. The legacy two_step option is rejected because it previously called svaseq.}
  \item{blas_threads}{Requested positive integer BLAS thread count. Actual threading depends on the BLAS library; this changes session thread settings. For SLURM, cannot exceed allocated CPUs.}
  \item{seed}{Random seed for reproducible feature diagnostics and SVA estimation.}
  \item{diagnostics}{Configuration argument; see the default in the usage and the module guide.}
  \item{diagnostics_dir}{Configuration argument; see the default in the usage and the module guide.}
  \item{return_selected_matrix}{Return the final centered/scaled M-value matrix used for fitting.}
  \item{return_feature_table}{Configuration argument; see the default in the usage and the module guide.}
  \item{return_logs}{Configuration argument; see the default in the usage and the module guide.}
  \item{write_rds}{Configuration argument; see the default in the usage and the module guide.}
  \item{rds_file}{Configuration argument; see the default in the usage and the module guide.}
  \item{log_file}{Configuration argument; see the default in the usage and the module guide.}
  \item{verbose}{Configuration argument; see the default in the usage and the module guide.}
  \item{.reuse_log}{Configuration argument; see the default in the usage and the module guide.}
}
\value{A preparation object, surrogate-variable result list, or SLURM job paths, depending on the function.}
\details{See the README and docs directory for input formats and optional dependencies.}
