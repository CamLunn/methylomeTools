# Default: Rscript inst/examples/install_dependencies.R
# Custom:  Rscript inst/examples/install_dependencies.R --lib /path/to/Rlibs
# With --lib, downloaded archives are retained in LIB/downloads.

install_methylome_dependencies <- function(lib = NULL) {
  install_args <- list()
  if (!is.null(lib)) {
    if (length(lib) != 1L || is.na(lib) || !nzchar(trimws(lib)))
      stop("lib must be a single nonempty directory path")
    lib <- path.expand(lib)
    if (!dir.exists(lib) && !dir.create(lib, recursive = TRUE, showWarnings = FALSE))
      stop("Could not create library directory: ", lib)
    lib <- normalizePath(lib, mustWork = TRUE)
    if (file.access(lib, 2L) != 0L) stop("Library directory is not writable: ", lib)
    downloads <- file.path(lib, "downloads")
    if (!dir.exists(downloads) &&
        !dir.create(downloads, recursive = TRUE, showWarnings = FALSE))
      stop("Could not create download directory: ", downloads)
    if (file.access(downloads, 2L) != 0L) stop("Download directory is not writable: ", downloads)
    .libPaths(c(lib, .libPaths()))
    install_args <- list(lib = lib, destdir = downloads)
    message("Installing new packages into: ", lib)
    message("Keeping downloaded archives in: ", downloads)
  } else {
    message("Using R's default library and download locations.")
  }

  repos <- getOption("repos")
  if (is.null(repos) || !length(repos)) repos <- c(CRAN = "https://cloud.r-project.org")
  if (!"CRAN" %in% names(repos)) repos <- c(repos, CRAN = "https://cloud.r-project.org")
  repos[repos == "@CRAN@"] <- "https://cloud.r-project.org"
  options(repos = repos)

  if (!requireNamespace("BiocManager", quietly = TRUE))
    do.call(utils::install.packages, c(list(pkgs = "BiocManager"), install_args))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    stop("BiocManager installation failed; inspect the installation output.")

  cran <- c("dplyr", "rlang", "ggplot2", "scales", "tibble", "tidyr", "testthat")
  bioc <- c("GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors", "BiocGenerics",
            "SummarizedExperiment", "bsseq", "sva", "DSS")
  do.call(BiocManager::install,
          c(list(pkgs = c(cran, bioc), ask = FALSE, update = FALSE), install_args))
  missing <- c(cran, bioc)[!vapply(c(cran, bioc), requireNamespace,
                                 logical(1), quietly = TRUE)]
  if (length(missing))
    stop("Required packages are still unavailable: ", paste(missing, collapse = ", "),
         ". Inspect the installation output.")
  message("Required packages are available. Install methylomeTools separately with R CMD INSTALL.")
  invisible(lib)
}

args <- commandArgs(trailingOnly = TRUE)
target_lib <- NULL
if (length(args)) {
  if (length(args) == 2L && identical(args[1L], "--lib")) {
    target_lib <- args[2L]
  } else if (length(args) == 1L && startsWith(args[1L], "--lib=")) {
    target_lib <- substring(args[1L], 7L)
  } else {
    stop("Usage: Rscript inst/examples/install_dependencies.R [--lib /path/to/Rlibs]")
  }
}
install_methylome_dependencies(lib = target_lib)

# Optional modules are installed separately. For custom libraries, set
# .libPaths(c('/path/to/Rlibs', .libPaths())) and pass lib='/path/to/Rlibs'
# and destdir='/path/to/Rlibs/downloads' to install.packages/BiocManager::install.
# BiocManager::install(c("ENmix", "patchwork", "mgcv", "ggrepel"))
# BiocManager::install(c("annotatr", "AnnotationDbi", "GenomicFeatures",
#   "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))
# install.packages(c("DBI", "RSQLite", "parallelly", "matrixStats", "RhpcBLASctl", "irlba"))
