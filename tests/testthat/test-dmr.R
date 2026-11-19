test_that("closed DMR boundaries map to the expected CpGs", {
  bs <- bsseq::BSseq(chr=rep("chr1", 4), pos=c(99,100,110,111),
    M=matrix(5L,4,3), Cov=matrix(10L,4,3), sampleNames=c("s1","s2","s3"))
  regions <- data.frame(chr="chr1", start_1based=100L, end_1based=110L)
  sites <- make_dmr_site_table(bs, regions)
  expect_identical(as.integer(sites$dms_pos), c(100L,110L))
  expect_equal(sites$dmr_width, c(11,11))
  long <- make_dmr_methylation_long(bs, sites, sites$dmr_id[1], min_coverage=1)
  expect_equal(unique(long$methylation_long$beta), .5)
})
