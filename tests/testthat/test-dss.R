test_that("DSS serial and PSOCK agree with direct chromosome fits", {
  skip_if_not_installed("DSS")
  set.seed(18)
  meta <- data.frame(sample_id=paste0("s", 1:16), exposure=rep(0:1, each=8))
  chr <- rep(c("chr1", "chr2"), each=100)
  pos <- rep(seq(100, by=20, length.out=100), 2)
  cov <- matrix(30L, 200, 16)
  meth <- matrix(rbinom(3200, 30, .4), 200, 16)
  bs <- bsseq::BSseq(chr=chr, pos=pos, M=meth, Cov=cov, sampleNames=meta$sample_id)
  run <- function(workers) {
    out <- tempfile("dss_")
    run_dss_parallel(bs, meta, ~ exposure, outcome_col="exposure", outcome_type="continuous",
      sample_id_col="sample_id", partition_mode="chromosome", n_workers=workers,
      parallel_backend="psock", outdir=out, result_prefix="test",
      save_compression="none", verbose=FALSE)
    tab <- readRDS(file.path(out, "test_merged_results.rds"))$main
    tab[order(tab$chr, tab$pos), ]
  }
  a <- run(1); b <- run(2)
  expect_equal(a$pvals, b$pvals, tolerance=1e-8)
  direct <- lapply(unique(chr), function(ch) {
    fit <- DSS::DMLfit.multiFactor(bs[chr == ch, ], design=meta, formula=~ exposure, smoothing=FALSE)
    DSS::DMLtest.multiFactor(fit, coef="exposure")
  })
  direct <- do.call(rbind, direct)
  direct <- direct[order(direct$chr, direct$pos), ]
  expect_equal(a$pvals, direct$pvals, tolerance=1e-8)
  expect_equal(a$fdr_global, p.adjust(a$pvals, "BH"))
})
