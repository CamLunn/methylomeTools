`%||%` <- function(x, y) if (is.null(x)) y else x

.require_optional <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("Install optional packages with BiocManager::install(): ",
                           paste(missing, collapse=", "), call.=FALSE)
  invisible(TRUE)
}

.available_workers <- function() {
  if (requireNamespace("parallelly", quietly=TRUE))
    return(max(1L, as.integer(parallelly::availableCores()) - 1L))
  allocated <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  if (is.na(allocated) || allocated < 1L) allocated <- 1L
  max(1L, allocated - 1L)
}
