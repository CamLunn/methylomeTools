test_that("SVA preparation accepts probe IDs and preserves direct sva results", {
  set.seed(11)
  x <- matrix(rnorm(500*24), 500, dimnames=list(paste0("cg", 1:500), paste0("s", 1:24)))
  meta <- data.frame(group=rep(0:1, each=12), row.names=colnames(x))
  mod <- model.matrix(~ group, meta); mod0 <- model.matrix(~ 1, meta)
  prep <- prepare_methyl_matrix(x, value_type="M", feature_strategy="hvf_only",
    n_hvf=200, override_hvf=TRUE, verbose=FALSE)
  result <- run_sva_on_matrix(prep, mod, mod0, num_sv_method="manual", num_sv_manual=1,
                              seed=3, verbose=FALSE)
  set.seed(3)
  direct <- sva::sva(prep$x_sva, mod, mod0, n.sv=1)
  expect_equal(abs(cor(result$sv$SV1, direct$sv[, 1])), 1, tolerance=1e-6)
  expect_identical(result$sv$sample_id, colnames(x))
  expect_error(run_sva_on_matrix(prep, mod[24:1, ], mod0,
    num_sv_method="manual", num_sv_manual=1, verbose=FALSE), "sample order")
})

test_that("beta preparation converts endpoints and preserves M-value inputs", {
  beta <- matrix(seq(0, 1, length.out = 30 * 12), 30,
                 dimnames = list(paste0("cg", 1:30), paste0("s", 1:12)))
  prepare <- function(x, type) {
    prepare_methyl_matrix(x, value_type = type, user_cpg_set = rownames(x),
      feature_strategy = "user_set_only", center_features = FALSE,
      scale_features = FALSE, verbose = FALSE)
  }
  converted <- beta_to_m(beta)
  from_beta <- prepare(beta, "beta")
  from_m <- prepare(converted, "M")
  expect_equal(from_beta$x_sva, converted)
  expect_equal(from_beta$x_sva, from_m$x_sva)
  expect_true(all(is.finite(from_beta$x_sva)))
})

test_that("SVA rejects invalid designs and stale beta preparations", {
  set.seed(5)
  x <- matrix(rnorm(40 * 12), 40,
              dimnames = list(paste0("cg", 1:40), paste0("s", 1:12)))
  meta <- data.frame(group = rep(0:1, each = 6), other = rep(0:2, 4),
                     row.names = colnames(x))
  mod <- model.matrix(~ group, meta)
  mod0 <- model.matrix(~ 1, meta)
  prep <- prepare_methyl_matrix(x, value_type = "M",
    user_cpg_set = rownames(x), feature_strategy = "user_set_only", verbose = FALSE)
  run <- function(full = mod, null = mod0, n = 1, object = prep) {
    run_sva_on_matrix(object, full, null, num_sv_method = "manual",
                      num_sv_manual = n, verbose = FALSE)
  }
  expect_error(run(cbind(mod, duplicate = mod[, 2])), "rank-deficient")
  expect_error(run(null = model.matrix(~ other, meta)), "nested")
  expect_error(run(n = 1.5), "non-negative integer")
  expect_error(run(n = 10), "allow at most")
  expect_error(run_sva_on_matrix(prep, mod, mod0, sva_method = "two_step"), "svaseq")
  prep$cfg$value_type <- "beta"
  prep$preparation_version <- NULL
  expect_error(run(object = prep), "Rebuild this cached beta preparation")
})

test_that("SVA worker requests describe actual process usage", {
  expect_warning(expect_identical(.sv_choose_workers(100, 1000, 4), 1L), "serial")
  expect_identical(.sv_choose_workers(100, 1000, NULL), 1L)
  expect_error(.sv_choose_workers(100, 1000, 1.5), "positive integer")
  expect_error(.sv_set_blas_threads(c(1, 2)), "positive integer")
})

test_that("SLURM SVA honors an explicit BLAS thread request", {
  out <- tempfile("sva_slurm_")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  x <- matrix(seq_len(120), 10,
              dimnames = list(paste0("cg", 1:10), paste0("s", 1:12)))
  meta <- data.frame(group = rep(0:1, each = 6), row.names = colnames(x))
  mod <- model.matrix(~ group, meta)
  mod0 <- model.matrix(~ 1, meta)
  submit_methyl_svs_slurm(x, mod, mod0, outdir = out, cpus_per_task = 8,
    blas_threads = 2, submit = FALSE, dry_run = TRUE)
  scripts <- list.files(out, pattern = "\\.(sh|sbatch)$", recursive = TRUE, full.names = TRUE)
  expect_gt(length(scripts), 0)
  script <- paste(unlist(lapply(scripts, readLines)), collapse = "\n")
  expect_match(script, "export OPENBLAS_NUM_THREADS=2", fixed = TRUE)
  rds <- list.files(out, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  objects <- lapply(rds, readRDS)
  expect_true(any(vapply(objects, function(x) is.list(x) && identical(x$blas_threads, 2L),
                         logical(1))))
})
