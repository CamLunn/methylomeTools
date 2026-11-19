test_that("beta conversion preserves scale, missingness and names", {
  b <- matrix(c(0, .25, .5, .75, 1, NA), 3, dimnames=list(letters[1:3], c("s1", "s2")))
  m <- beta_to_m(b)
  expect_equal(m[2, 1], log2(.25/.75))
  expect_true(is.finite(m[1, 1]))
  expect_true(is.na(m[3, 2]))
  expect_identical(dimnames(m), dimnames(b))
  expect_error(beta_to_m(matrix(1.2)), "between")
})

test_that("array mappings retain row order and handle absent probes", {
  b <- matrix(c(.1, .3, .5, .2, .4, .6), 3,
              dimnames=list(c("p1", "p2", "p3"), c("s1", "s2")))
  map <- data.frame(cpg_id=c("cgA", "cgB", "cgC"), probe_ids=c("p3,p1,p2", "absent", "p2"))
  x <- align_array_beta_to_cpgs(b, map)
  expect_identical(rownames(x), map$cpg_id)
  expect_equal(as.numeric(x[1, ]), c(.3, .4))
  expect_true(all(is.na(x[2, ])))
  expect_equal(as.numeric(x[3, ]), c(.3, .4))
})
