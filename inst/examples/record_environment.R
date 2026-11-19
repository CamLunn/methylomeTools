outdir <- "outputs/environment"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
if (requireNamespace("BiocManager", quietly=TRUE))
  writeLines(as.character(BiocManager::version()), file.path(outdir, "bioconductor_version.txt"))
pkgs <- installed.packages()
utils::write.csv(pkgs[, c("Package", "Version", "LibPath")],
                 file.path(outdir, "installed_packages.csv"), row.names=FALSE)
