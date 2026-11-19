# Methylation feature preparation, SVA fitting and SLURM submission.

build_methyl_svs <- function(
    meth_mat,
    mod,
    mod0,
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
    num_sv_method = c("be", "leek", "manual"),
    num_sv_manual = NULL,
    sva_method = c("irw", "two_step"),
    center_features = TRUE,
    scale_features = FALSE,
    use_delayed = TRUE,
    workers = NULL,
    blas_threads = NULL,
    seed = 1,
    log_file = NULL,
    verbose = TRUE,
    diagnostics = NULL,
    diagnostics_dir = NULL,
    return_selected_matrix = FALSE,
    return_feature_table = TRUE,
    return_logs = TRUE,
    write_rds = FALSE,
    rds_file = NULL,
    cache_prep_rds = NULL
) {
  value_type       <- match.arg(value_type)
  feature_strategy <- match.arg(feature_strategy)
  hvf_stat         <- match.arg(hvf_stat)
  impute_method    <- match.arg(impute_method)
  num_sv_method    <- match.arg(num_sv_method)
  sva_method       <- match.arg(sva_method)
  if (sva_method == "two_step") {
    stop("The legacy two_step option called svaseq. Use sva_method = 'irw' for methylation M values.")
  }
  
  prep <- prepare_methyl_matrix(
    meth_mat = meth_mat,
    value_type = value_type,
    user_cpg_set = user_cpg_set,
    feature_strategy = feature_strategy,
    n_hvf = n_hvf,
    override_hvf = override_hvf,
    hvf_stat = hvf_stat,
    missingness_max = missingness_max,
    impute = impute,
    impute_method = impute_method,
    knn_k = knn_k,
    center_features = center_features,
    scale_features = scale_features,
    use_delayed = use_delayed,
    workers = workers,
    blas_threads = blas_threads,
    seed = seed,
    log_file = log_file,
    verbose = verbose,
    write_rds = !is.null(cache_prep_rds),
    rds_file = cache_prep_rds,
    .keep_log_open = TRUE          # let run stage append to the same log
  )
  
  out <- run_sva_on_matrix(
    prep = prep,
    mod = mod,
    mod0 = mod0,
    num_sv_method = num_sv_method,
    num_sv_manual = num_sv_manual,
    sva_method = sva_method,
    blas_threads = blas_threads,
    seed = seed,
    diagnostics = diagnostics,
    diagnostics_dir = diagnostics_dir,
    return_selected_matrix = return_selected_matrix,
    return_feature_table = return_feature_table,
    return_logs = return_logs,
    write_rds = write_rds,
    rds_file = rds_file,
    verbose = verbose,
    .reuse_log = TRUE
  )
  
  out
}


#### STAGE 1: prepare_methyl_matrix()

prepare_methyl_matrix <- function(
    meth_mat,
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
    .keep_log_open = FALSE
) {
  value_type       <- match.arg(value_type)
  feature_strategy <- match.arg(feature_strategy)
  hvf_stat         <- match.arg(hvf_stat)
  impute_method    <- match.arg(impute_method)
  
  .svlog_init(log_file = log_file, verbose = verbose)
  if (!isTRUE(.keep_log_open)) on.exit(.svlog_close(), add = TRUE)
  
  t0 <- Sys.time()
  set.seed(seed)
  
  .sv_set_blas_threads(blas_threads, verbose = verbose)
  
  .svlog("Starting prepare_methyl_matrix()")
  .svlog(sprintf("Input class: %s", paste(class(meth_mat), collapse = ", ")))
  .svlog(sprintf("Input dims: %d CpGs x %d samples", nrow(meth_mat), ncol(meth_mat)))
  .svlog(sprintf("value_type=%s, feature_strategy=%s, hvf_stat=%s",
                 value_type, feature_strategy, hvf_stat))
  
  cfg <- .sv_build_config(
    value_type = value_type,
    feature_strategy = feature_strategy,
    n_hvf = n_hvf,
    override_hvf = override_hvf,
    hvf_stat = hvf_stat,
    missingness_max = missingness_max,
    impute = impute,
    impute_method = impute_method,
    knn_k = knn_k,
    center_features = center_features,
    scale_features = scale_features,
    use_delayed = use_delayed,
    workers = workers,
    blas_threads = blas_threads
  )
  
  .sv_validate_prep_inputs(
    meth_mat = meth_mat,
    user_cpg_set = user_cpg_set,
    cfg = cfg
  )
  
  cfg$workers <- .sv_choose_workers(
    n_samples = ncol(meth_mat),
    n_cpg = nrow(meth_mat),
    workers = cfg$workers
  )
  cfg$requested_workers <- workers
  .svlog("Preprocessing runs in one R process; workers does not create a pool.")
  
  .svlog("Parsing CpG rownames")
  cpg_info <- .sv_time_it(
    "CpG parsing",
    .sv_parse_cpg_ids(rownames(meth_mat))
  )
  
  .svlog("Preparing user CpG set")
  user_set_info <- .sv_time_it(
    "User CpG preparation",
    .sv_prepare_user_set(
      all_cpg_ids = rownames(meth_mat),
      user_cpg_set = user_cpg_set
    )
  )
  .svlog(sprintf(
    "User CpGs supplied=%d, present_in_matrix=%d",
    user_set_info$n_user_supplied,
    user_set_info$n_user_present
  ))
  
  .svlog("Computing missingness statistics")
  missing_stats <- .sv_time_it(
    "Missingness computation",
    .sv_compute_missingness(
      meth_mat = meth_mat,
      use_delayed = cfg$use_delayed,
      workers = cfg$workers
    )
  )
  
  .svlog("Selecting feature set for SVA")
  feature_plan <- .sv_time_it(
    "Feature selection",
    .sv_select_features(
      meth_mat = meth_mat,
      cpg_info = cpg_info,
      user_set_info = user_set_info,
      missing_stats = missing_stats,
      feature_strategy = cfg$feature_strategy,
      n_hvf = cfg$n_hvf,
      override_hvf = cfg$override_hvf,
      hvf_stat = cfg$hvf_stat,
      missingness_max = cfg$missingness_max,
      use_delayed = cfg$use_delayed,
      workers = cfg$workers
    )
  )
  
  .svlog(sprintf(
    "Selected %d CpGs total (%d user-set, %d HVF)",
    feature_plan$summary$n_selected_total,
    feature_plan$summary$n_user_set_used,
    feature_plan$summary$n_hvf_used
  ))
  
  if (feature_plan$summary$n_selected_total < 2L) {
    stop("Too few CpGs selected for SVA after filtering.")
  }
  
  .svlog("Extracting selected matrix into memory")
  x_sel <- .sv_time_it(
    "Matrix extraction",
    .sv_extract_selected_matrix(
      meth_mat = meth_mat,
      selected_cpg_ids = feature_plan$selected_cpg_ids
    )
  )
  
  if (cfg$impute) {
    .svlog("Imputation requested")
    x_sel <- .sv_time_it(
      "Imputation",
      .sv_impute_selected_matrix(
        x = x_sel,
        cpg_info = cpg_info[match(rownames(x_sel), cpg_info$cpg_id), , drop = FALSE],
        method = cfg$impute_method,
        knn_k = cfg$knn_k,
        workers = cfg$workers
      )
    )
  } else {
    .svlog("Imputation disabled; removing any rows with remaining missing values")
    keep_complete <- rowSums(is.na(x_sel)) == 0
    n_drop <- sum(!keep_complete)
    if (n_drop > 0) {
      .svlog(sprintf("Dropping %d selected CpGs with remaining NA values", n_drop))
      x_sel <- x_sel[keep_complete, , drop = FALSE]
    }
  }
  
  if (nrow(x_sel) < 2L) {
    stop("Too few complete CpGs remain after NA handling to run SVA.")
  }
  
  if (cfg$value_type == "beta") {
    .svlog("Converting selected beta values to log2 M values (eps = 1e-6)")
    x_sel <- beta_to_m(x_sel)
  }
  cfg$matrix_scale <- "M"

  .svlog("Preparing selected matrix for SVA")
  x_sva <- .sv_time_it(
    "Matrix preparation for SVA",
    .sv_prepare_matrix_for_sva(
      x = x_sel,
      center_features = cfg$center_features,
      scale_features = cfg$scale_features
    )
  )
  
  prep_runtime <- list(
    start_time = as.character(t0),
    end_time = as.character(Sys.time()),
    elapsed_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    input_dim = c(nrow(meth_mat), ncol(meth_mat)),
    selected_dim = c(nrow(x_sel), ncol(x_sel)),
    sva_dim = c(nrow(x_sva), ncol(x_sva)),
    workers = cfg$workers
  )
  
  prep <- list(
    x_sva = x_sva,
    x_sel = x_sel,
    selected_cpgs = rownames(x_sel),
    sample_ids = colnames(meth_mat),
    preparation_version = 2L,
    feature_plan = feature_plan,
    cpg_info = cpg_info,
    missing_stats = missing_stats,
    cfg = cfg,
    prep_runtime = prep_runtime,
    logs = .svlog_get()
  )
  class(prep) <- c("methyl_prep", "list")
  
  if (isTRUE(write_rds)) {
    if (is.null(rds_file)) stop("write_rds=TRUE but rds_file is NULL")
    .svlog(sprintf("Caching prepared matrix to: %s", rds_file))
    dir.create(dirname(rds_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(prep, rds_file)
  }
  
  .svlog("Finished prepare_methyl_matrix()")
  prep
}


#### STAGE 2: run_sva_on_matrix()

run_sva_on_matrix <- function(
    prep,
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
    .reuse_log = FALSE
) {
  num_sv_method <- match.arg(num_sv_method)
  sva_method    <- match.arg(sva_method)
  if (sva_method == "two_step") {
    stop("The legacy two_step option called svaseq, a sequencing-data estimator. ",
         "Use sva_method = 'irw' for methylation M values.")
  }
  
  if (is.character(prep) && length(prep) == 1L) {
    prep <- readRDS(prep)
  }
  if (!inherits(prep, "methyl_prep")) {
    stop("prep must be a methyl_prep object (from prepare_methyl_matrix) or a path to one.")
  }
  
  if (!isTRUE(.reuse_log)) {
    .svlog_init(log_file = log_file, verbose = verbose)
    on.exit(.svlog_close(), add = TRUE)
    .svlog("Starting run_sva_on_matrix()")
  } else {
    on.exit(.svlog_close(), add = TRUE)
  }
  
  t0 <- Sys.time()
  set.seed(seed)
  .sv_set_blas_threads(blas_threads, verbose = verbose)
  
  x_sva         <- prep$x_sva
  x_sel         <- prep$x_sel
  feature_plan  <- prep$feature_plan
  cpg_info      <- prep$cpg_info
  missing_stats <- prep$missing_stats
  sample_names  <- prep$sample_ids
  cfg           <- prep$cfg
  if (identical(cfg$value_type, "beta") &&
      !identical(prep$preparation_version, 2L)) {
    stop("Rebuild this cached beta preparation with prepare_methyl_matrix(); ",
         "older versions did not convert beta values to M values.")
  }
  cfg$workers <- 1L
  cfg$blas_threads <- blas_threads
  
  cfg$num_sv_method          <- num_sv_method
  cfg$num_sv_manual          <- num_sv_manual
  cfg$sva_method             <- sva_method
  cfg$diagnostics            <- diagnostics
  cfg$return_selected_matrix <- return_selected_matrix
  cfg$return_feature_table   <- return_feature_table
  cfg$return_logs            <- return_logs
  cfg$write_rds              <- write_rds
  cfg$rds_file               <- rds_file
  
  for (design in list(mod, mod0)) {
    if (!is.null(rownames(design)) &&
        !identical(rownames(design), as.character(seq_len(nrow(design)))) &&
        !identical(rownames(design), colnames(prep$x_sva)))
      stop("Named model-matrix rows must match methylation sample order")
  }
  .sv_validate_model_inputs(mod = mod, mod0 = mod0, n_samples = length(sample_names), cfg = cfg)
  
  if (!is.null(diagnostics) && length(diagnostics) > 0) {
    diagnostics_dir <- .sv_prepare_diagnostics_dir(diagnostics_dir)
    cfg$diagnostics_dir <- diagnostics_dir
    .svlog(sprintf("Diagnostics will be written to: %s", diagnostics_dir))
  }
  
  .svlog("Estimating number of surrogate variables")
  n_sv <- .sv_time_it(
    "num.sv estimation",
    .sv_choose_num_sv(
      x = x_sva,
      mod = mod,
      mod0 = mod0,
      method = cfg$num_sv_method,
      manual_n = cfg$num_sv_manual
    )
  )
  max_sv <- max(0L, min(nrow(x_sva), ncol(x_sva) - qr(mod)$rank - 1L))
  if (n_sv > max_sv) {
    stop("Requested ", n_sv, " SVs, but this matrix and design allow at most ", max_sv,
         " while retaining residual degrees of freedom.")
  }
  .svlog(sprintf("n_sv selected: %d", n_sv))
  
  if (n_sv < 1L) {
    warning("Estimated n_sv < 1; returning zero SVs.")
    sv_df <- data.frame(sample_id = sample_names, stringsAsFactors = FALSE)
    diagnostics_out <- NULL
  } else {
    .svlog("Running SVA")
    sva_fit <- .sv_time_it(
      "SVA fit",
      .sv_run_sva(
        x = x_sva,
        mod = mod,
        mod0 = mod0,
        n_sv = n_sv,
        method = cfg$sva_method
      )
    )
    
    .svlog("Formatting SVA output")
    sv_df <- .sv_format_sv_output(
      sva_fit = sva_fit,
      sample_names = sample_names
    )
    
    diagnostics_out <- NULL
    if (!is.null(cfg$diagnostics) && length(cfg$diagnostics) > 0) {
      .svlog("Generating requested diagnostics")
      diagnostics_out <- .sv_time_it(
        "Diagnostics generation",
        .sv_generate_diagnostics(
          x = x_sva,
          sv_df = sv_df,
          feature_plan = feature_plan,
          cpg_info = cpg_info,
          missing_stats = missing_stats,
          diagnostics = cfg$diagnostics,
          diagnostics_dir = cfg$diagnostics_dir
        )
      )
    }
  }
  
  runtime <- .sv_collect_runtime_summary(
    start_time = t0,
    input_dim = prep$prep_runtime$input_dim,
    x_sel = x_sel,
    x_sva = x_sva,
    feature_plan = feature_plan,
    n_sv = n_sv,
    workers = cfg$workers
  )
  runtime$prep_elapsed_seconds <- prep$prep_runtime$elapsed_seconds
  
  out <- list(
    sv = sv_df,
    n_sv = n_sv,
    selected_cpgs = rownames(x_sel),
    feature_summary = feature_plan$summary,
    settings = cfg,
    runtime = runtime,
    diagnostics = diagnostics_out
  )
  
  if (cfg$return_feature_table) {
    ft <- feature_plan$feature_table
    ft$kept_after_final_na_handling <- ft$cpg_id %in% rownames(x_sel)
    out$feature_table <- ft
  }
  
  if (cfg$return_selected_matrix) {
    out$selected_matrix <- x_sva
  }
  
  if (cfg$return_logs) {
    out$logs <- .svlog_get()
  }
  
  if (isTRUE(cfg$write_rds)) {
    if (is.null(cfg$rds_file)) stop("write_rds=TRUE but rds_file is NULL")
    .svlog(sprintf("Writing output RDS to: %s", cfg$rds_file))
    dir.create(dirname(cfg$rds_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cfg$rds_file)
  }
  
  .svlog("Finished run_sva_on_matrix()")
  out
}


#### SLURM launcher

submit_methyl_svs_slurm <- function(
    meth_mat,
    mod,
    mod0,
    outdir,
    ...,                                   # passed straight to build_methyl_svs()
    sva_lib = NULL,
    staging_dir = file.path(outdir, "slurm_staging"),
    job_name = "methyl_sva",
    cpus_per_task = 16L,
    mem_gb = NULL,                         # auto-estimated if NULL
    time = "12:00:00",
    partition = NULL,                      # e.g. "fn_medium"; NULL -> cluster default
    r_module = NULL,
    extra_modules = character(0),
    openblas_so = NULL,                    # optional LD_PRELOAD path to a threaded BLAS
    blas_threads = NULL,                   # defaults to cpus_per_task on the node
    account = NULL,
    email = NULL,
    sbatch_bin = NULL,
    submit = TRUE,
    dry_run = FALSE
) {
  .sv_positive_integer(cpus_per_task, "cpus_per_task")
  if (is.null(sva_lib)) sva_lib <- "-"
  stopifnot(is.character(sva_lib), length(sva_lib) == 1L)
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
  
  build_args <- list(...)
  # The SVA estimator runs in one process; BLAS may use several threads.
  if (is.null(build_args$workers))      build_args$workers      <- 1L
  if (is.null(build_args$use_delayed))  build_args$use_delayed  <- FALSE
  if (!is.null(blas_threads)) {
    .sv_positive_integer(blas_threads, "blas_threads")
    if (blas_threads > cpus_per_task) stop("blas_threads cannot exceed cpus_per_task")
  }
  build_args$blas_threads <- if (is.null(blas_threads)) {
    "${SLURM_CPUS_PER_TASK}"
  } else {
    as.integer(blas_threads)
  }
  
  # Save batch results unless the caller explicitly opts out.
  if (is.null(build_args$write_rds)) build_args$write_rds <- TRUE
  if (isTRUE(build_args$write_rds) && is.null(build_args$rds_file)) {
    build_args$rds_file <- file.path(outdir, "build_methyl_svs_result.rds")
  }
  # Cache preparation for models using the same samples and settings.
  if (is.null(build_args$cache_prep_rds)) {
    build_args$cache_prep_rds <- file.path(outdir, "prep.rds")
  }
  result_rds <- build_args$rds_file
  if (isFALSE(build_args$write_rds)) {
    warning("write_rds = FALSE: the SLURM job will print n_sv but save NO result file. ",
            "You will not be able to read the result back. Set write_rds = TRUE (default).")
    result_rds <- NA_character_
  }
  
  res <- .sv_estimate_slurm_resources(meth_mat = meth_mat, mem_gb = mem_gb)
  mem_gb <- res$mem_gb
  
  inputs_path <- file.path(staging_dir, "inputs.rds")
  args_path   <- file.path(staging_dir, "build_args.rds")
  driver_path <- file.path(staging_dir, "_sva_driver.R")
  sbatch_path <- file.path(staging_dir, "_sva_job.sbatch")
  
  message(sprintf("Serializing inputs to %s (this can take a while for large matrices)...", inputs_path))
  saveRDS(list(meth_mat = meth_mat, mod = mod, mod0 = mod0), inputs_path)
  saveRDS(build_args, args_path)
  
  .sv_write_slurm_driver(driver_path)
  .sv_write_sbatch(
    sbatch_path   = sbatch_path,
    driver_path   = driver_path,
    sva_lib       = sva_lib,
    inputs_path   = inputs_path,
    args_path     = args_path,
    outdir        = outdir,
    job_name      = job_name,
    cpus_per_task = cpus_per_task,
    mem_gb        = mem_gb,
    time          = time,
    partition     = partition,
    r_module      = r_module,
    extra_modules = extra_modules,
    openblas_so   = openblas_so,
    blas_threads  = blas_threads,
    account       = account,
    email         = email
  )
  
  message(sprintf("Wrote sbatch script: %s", sbatch_path))
  message(sprintf("Estimated resources: %d CPUs, %d GB, time %s", cpus_per_task, mem_gb, time))
  if (!is.na(result_rds)) {
    message(sprintf("On completion, read results with:  readRDS(\"%s\")", result_rds))
  }
  
  job_id <- NA_character_
  if (isTRUE(submit) && !isTRUE(dry_run)) {
    sbatch_resolved <- .sv_find_sbatch(sbatch_bin)
    if (!is.null(sbatch_resolved)) {
      out <- suppressWarnings(system2(sbatch_resolved, sbatch_path, stdout = TRUE, stderr = TRUE))
      message(paste(out, collapse = "\n"))
      m <- regmatches(out, regexpr("[0-9]+", out))
      if (length(m)) job_id <- m[[1]]
    } else {
      warning("Could not locate sbatch; not submitting. Submit manually:\n  sbatch ", sbatch_path,
              "\nOr pass sbatch_bin='/path/to/sbatch' to submit_methyl_svs_slurm().")
    }
  } else {
    message("dry_run/submit=FALSE: not submitting. To submit:\n  sbatch ", sbatch_path)
  }
  
  invisible(list(
    job_id = job_id,
    result_rds = result_rds,
    prep_rds = build_args$cache_prep_rds,
    sbatch = sbatch_path,
    driver = driver_path,
    inputs = inputs_path,
    build_args = args_path,
    outdir = outdir,
    mem_gb = mem_gb,
    cpus_per_task = cpus_per_task
  ))
}

# Locate sbatch: explicit arg -> PATH -> known cluster locations.
.sv_find_sbatch <- function(sbatch_bin = NULL) {
  if (!is.null(sbatch_bin) && nzchar(sbatch_bin) && file.exists(sbatch_bin)) {
    return(sbatch_bin)
  }
  found <- Sys.which("sbatch")
  if (nzchar(found)) return(unname(found))
  
  candidates <- c(
    "/cm/shared/apps/slurm/current/bin/sbatch",  # BigPurple
    "/usr/bin/sbatch",
    "/usr/local/bin/sbatch",
    "/opt/slurm/bin/sbatch"
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(hit[[1]])
  
  NULL
}

.sv_estimate_slurm_resources <- function(meth_mat, mem_gb = NULL) {
  input_gb <- (as.numeric(nrow(meth_mat)) * as.numeric(ncol(meth_mat)) * 8) / 1e9
  if (is.null(mem_gb)) {
    # raw input + filtered copy + working copies during centering/SVA.
    mem_gb <- max(64L, as.integer(ceiling(input_gb * 6)))
  }
  list(input_gb = input_gb, mem_gb = mem_gb)
}

.sv_write_slurm_driver <- function(driver_path) {
  driver <- '#!/usr/bin/env Rscript
## Auto-generated by submit_methyl_svs_slurm(). Runs on the compute node.
args <- commandArgs(trailingOnly = TRUE)
sva_lib     <- args[[1]]
inputs_path <- args[[2]]
args_path   <- args[[3]]

## Resolve BLAS threads from the SLURM allocation.
barg <- readRDS(args_path)
threads <- if (identical(barg$blas_threads, "${SLURM_CPUS_PER_TASK}") ||
               is.null(barg$blas_threads)) {
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
} else {
  as.integer(barg$blas_threads)
}
Sys.setenv(OMP_NUM_THREADS = threads,
           OPENBLAS_NUM_THREADS = threads,
           MKL_NUM_THREADS = threads)

if (identical(sva_lib, "-")) library(methylomeTools) else source(sva_lib)

inp  <- readRDS(inputs_path)
## Replace the threads placeholder with the real allocation.
if (identical(barg$blas_threads, "${SLURM_CPUS_PER_TASK}") || is.null(barg$blas_threads)) {
  barg$blas_threads <- threads
}

call_args <- c(list(meth_mat = inp$meth_mat, mod = inp$mod, mod0 = inp$mod0), barg)
res <- do.call(build_methyl_svs, call_args)
cat(sprintf("Done. n_sv = %s\\n", res$n_sv))
'
writeLines(driver, driver_path)
Sys.chmod(driver_path, mode = "0755")
invisible(driver_path)
}

.sv_write_sbatch <- function(
    sbatch_path, driver_path, sva_lib, inputs_path, args_path, outdir,
    job_name, cpus_per_task, mem_gb, time, partition,
    r_module, extra_modules, openblas_so, blas_threads, account, email
) {
  L <- c()
  add <- function(...) L <<- c(L, sprintf(...))
  
  add("#!/bin/bash")
  add("#SBATCH --job-name=%s", job_name)
  add("#SBATCH --nodes=1")
  add("#SBATCH --ntasks=1")
  add("#SBATCH --cpus-per-task=%d", as.integer(cpus_per_task))
  add("#SBATCH --mem=%dG", as.integer(mem_gb))
  add("#SBATCH --time=%s", time)
  add("#SBATCH --output=%s", file.path(outdir, "logs", "sva_%j.out"))
  add("#SBATCH --error=%s",  file.path(outdir, "logs", "sva_%j.err"))


  if (!is.null(partition)) add("#SBATCH --partition=%s", partition)
  if (!is.null(account))   add("#SBATCH --account=%s", account)
  if (!is.null(email)) {
    add("#SBATCH --mail-type=END,FAIL")
    add("#SBATCH --mail-user=%s", email)
  }
  add("")
  add("set -eo pipefail")   # not -u: module() often references unbound vars
  add("")
  add("# `module` is a shell function defined by an init script that login shells")
  add("# source but a non-interactive sbatch shell does not. Source it first, or")
  add("# `module purge` fails with 'command not found' and -e aborts the job.")
  add("if ! command -v module >/dev/null 2>&1; then")
  add("  for f in /etc/profile.d/modules.sh \\")
  add("           /usr/share/Modules/init/bash \\")
  add("           /cm/local/apps/environment-modules/current/init/bash; do")
  add("    [ -r \"$f\" ] && source \"$f\" && break")
  add("  done")
  add("fi")
  add("")
  if (!is.null(r_module)) {
    add("module purge")
    add("module load %s", shQuote(r_module))
  }
  for (m in extra_modules) add("module load %s", shQuote(m))
  add("")
  thread_value <- if (is.null(blas_threads)) "$SLURM_CPUS_PER_TASK" else as.character(blas_threads)
  add("# Set threads before R starts; actual threading depends on the BLAS library.")
  add("export OMP_NUM_THREADS=%s", thread_value)
  add("export OPENBLAS_NUM_THREADS=%s", thread_value)
  add("export MKL_NUM_THREADS=%s", thread_value)
  if (!is.null(openblas_so)) {
    add("")
    add("# Reference BLAS is single-threaded. Preload a threaded BLAS so the")
    add("# SVDs inside num.sv / sva use all allocated cores.")
    add("export LD_PRELOAD=%s", shQuote(openblas_so))
  } else if (!is.null(r_module)) {
    add("")
    add("# NOTE: if `module load %s` gives you reference BLAS, the cores above", r_module)
    add("# will NOT be used. Either load an R module built against MKL/OpenBLAS,")
    add("# or set openblas_so= in submit_methyl_svs_slurm() to LD_PRELOAD one.")
  }
  add("")
  add("Rscript %s \\", shQuote(driver_path))
  add("  %s \\", shQuote(sva_lib))
  add("  %s \\", shQuote(inputs_path))
  add("  %s", shQuote(args_path))
  add("")
  
  writeLines(L, sbatch_path)
  invisible(sbatch_path)
}


#### BLAS thread control

.sv_blas_info <- function() {
  blas <- tryCatch(extSoftVersion()[["BLAS"]], error = function(e) NA_character_)
  is_ref <- is.na(blas) || grepl("libRblas", blas, fixed = TRUE)
  list(path = blas, is_reference = is_ref)
}

.sv_set_blas_threads <- function(n_threads = NULL, verbose = TRUE) {
  if (is.null(n_threads)) return(invisible(NULL))
  .sv_positive_integer(n_threads, "blas_threads")
  n_threads <- as.integer(n_threads)
  
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(n_threads), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(n_threads), silent = TRUE)
  }
  if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    warning("RhpcBLASctl is unavailable; thread environment variables may not ",
            "reconfigure an already loaded BLAS. Set them before starting R.")
  }
  Sys.setenv(
    OMP_NUM_THREADS = n_threads,
    OPENBLAS_NUM_THREADS = n_threads,
    MKL_NUM_THREADS = n_threads
  )
  
  info <- .sv_blas_info()
  if (isTRUE(verbose)) {
    if (info$is_reference) {
      .svlog(sprintf(
        "WARNING: BLAS looks like reference single-threaded (%s). Requested %d threads may not be available; check the BLAS configuration.",
        info$path, n_threads))
    } else {
      .svlog(sprintf("BLAS=%s; requested %d threads (actual use depends on the library)", info$path, n_threads))
    }
  }
  invisible(NULL)
}


#### Logging helpers

.svlog_env <- new.env(parent = emptyenv())

.svlog_init <- function(log_file = NULL, verbose = TRUE) {
  .svlog_env$messages <- character()
  .svlog_env$log_file <- log_file
  .svlog_env$verbose <- isTRUE(verbose)
  
  if (!is.null(log_file)) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    cat("", file = log_file, append = FALSE)
  }
}

.svlog <- function(msg) {
  if (is.null(.svlog_env$messages)) .svlog_env$messages <- character()
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", stamp, msg)
  .svlog_env$messages <- c(.svlog_env$messages, line)
  
  if (isTRUE(.svlog_env$verbose)) {
    message(line)
  }
  if (!is.null(.svlog_env$log_file)) {
    cat(line, "\n", file = .svlog_env$log_file, append = TRUE)
  }
}

.svlog_get <- function() {
  .svlog_env$messages
}

.svlog_close <- function() {
  invisible(NULL)
}

.sv_time_it <- function(label, expr) {
  t0 <- Sys.time()
  out <- expr
  t1 <- Sys.time()
  .svlog(sprintf("%s completed in %.2f sec", label, as.numeric(difftime(t1, t0, units = "secs"))))
  out
}


#### Config + validation

.sv_build_config <- function(
    value_type,
    feature_strategy,
    n_hvf,
    override_hvf,
    hvf_stat,
    missingness_max,
    impute,
    impute_method,
    knn_k,
    center_features,
    scale_features,
    use_delayed,
    workers,
    blas_threads
) {
  list(
    value_type = value_type,
    feature_strategy = feature_strategy,
    n_hvf = n_hvf,
    override_hvf = override_hvf,
    hvf_stat = hvf_stat,
    missingness_max = missingness_max,
    impute = impute,
    impute_method = impute_method,
    knn_k = knn_k,
    center_features = center_features,
    scale_features = scale_features,
    use_delayed = use_delayed,
    workers = workers,
    blas_threads = blas_threads
  )
}

.sv_validate_prep_inputs <- function(meth_mat, user_cpg_set, cfg) {
  if (isTRUE(cfg$impute)) stop("Imputation is not implemented; use impute = FALSE")
  if (is.matrix(meth_mat) && !is.numeric(meth_mat)) stop("meth_mat must be numeric")
  if (!(is.matrix(meth_mat) || methods::is(meth_mat, "DelayedMatrix"))) {
    stop("meth_mat must be a matrix or a DelayedMatrix")
  }
  if (nrow(meth_mat) < 2L || ncol(meth_mat) < 2L) {
    stop("meth_mat must have at least 2 rows and 2 columns")
  }
  if (is.null(rownames(meth_mat))) {
    stop("meth_mat must have unique feature IDs as rownames")
  }
  if (is.null(colnames(meth_mat))) {
    stop("meth_mat must have sample IDs as colnames")
  }
  if (anyDuplicated(rownames(meth_mat))) {
    stop("meth_mat rownames must be unique")
  }
  if (anyDuplicated(colnames(meth_mat))) {
    stop("meth_mat colnames must be unique")
  }
  
  if (!is.null(user_cpg_set) && !is.character(user_cpg_set)) {
    stop("user_cpg_set must be NULL or a character vector")
  }
  
  if (!is.numeric(cfg$missingness_max) || length(cfg$missingness_max) != 1L ||
      is.na(cfg$missingness_max) || cfg$missingness_max < 0 || cfg$missingness_max > 1) {
    stop("missingness_max must be a single numeric value in [0, 1]")
  }
  if (!is.numeric(cfg$knn_k) || length(cfg$knn_k) != 1L || is.na(cfg$knn_k) || cfg$knn_k < 1) {
    stop("knn_k must be a positive integer")
  }
  if (!is.null(cfg$n_hvf)) {
    if (!is.numeric(cfg$n_hvf) || length(cfg$n_hvf) != 1L || is.na(cfg$n_hvf) || cfg$n_hvf < 1) {
      stop("n_hvf must be NULL or a positive integer")
    }
  }
  invisible(TRUE)
}

.sv_validate_model_inputs <- function(mod, mod0, n_samples, cfg) {
  if (!is.matrix(mod) || !is.matrix(mod0)) {
    stop("mod and mod0 must be model matrices")
  }
  if (nrow(mod) != n_samples) {
    stop("nrow(mod) must equal the number of samples in the prepared matrix")
  }
  if (nrow(mod0) != n_samples) {
    stop("nrow(mod0) must equal the number of samples in the prepared matrix")
  }
  if (nrow(mod) != nrow(mod0)) {
    stop("mod and mod0 must have same number of rows")
  }
  if (!is.numeric(mod) || !is.numeric(mod0)) stop("Model matrices must be numeric")
  if (anyNA(mod) || anyNA(mod0) || any(!is.finite(mod)) || any(!is.finite(mod0)))
    stop("Model matrices must contain finite values only")
  if (ncol(mod) < 1L || ncol(mod0) < 1L) stop("Model matrices must have columns")
  full_qr <- qr(mod)
  if (full_qr$rank < ncol(mod)) stop("mod is rank-deficient")
  if (qr(mod0)$rank < ncol(mod0)) stop("mod0 is rank-deficient")
  if (ncol(mod0) >= ncol(mod)) stop("mod0 must be a smaller model nested within mod")
  if (n_samples - full_qr$rank < 2L) stop("mod leaves too few residual degrees of freedom")
  # Compare column spaces, allowing equivalent encodings of the null model.
  residual <- qr.resid(full_qr, mod0)
  if (sqrt(sum(residual^2)) > 1e-7 * max(1, sqrt(sum(mod0^2)))) {
    stop("mod0 must be nested within mod")
  }
  
  if (cfg$num_sv_method == "manual") {
    if (is.null(cfg$num_sv_manual)) {
      stop("num_sv_manual must be supplied when num_sv_method='manual'")
    }
    if (!is.numeric(cfg$num_sv_manual) || length(cfg$num_sv_manual) != 1L ||
        !is.finite(cfg$num_sv_manual) || cfg$num_sv_manual < 0 ||
        cfg$num_sv_manual != floor(cfg$num_sv_manual) ||
        cfg$num_sv_manual > .Machine$integer.max) {
      stop("num_sv_manual must be a non-negative integer")
    }
  }
  
  if (!is.null(cfg$diagnostics)) {
    valid_diags <- c("feature_summary", "missingness", "sv_cor", "scree", "sv_pairs", "runtime")
    bad <- setdiff(cfg$diagnostics, valid_diags)
    if (length(bad) > 0) {
      stop(sprintf("Unsupported diagnostics requested: %s", paste(bad, collapse = ", ")))
    }
  }
  
  invisible(TRUE)
}

.sv_prepare_diagnostics_dir <- function(diagnostics_dir = NULL) {
  if (is.null(diagnostics_dir)) {
    diagnostics_dir <- file.path(getwd(), "sva_diagnostics")
  }
  dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
  diagnostics_dir
}


#### Worker logic

.sv_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop(name, " must be a positive integer")
  }
  invisible(TRUE)
}

.sv_choose_workers <- function(n_samples, n_cpg, workers = NULL) {
  if (!is.null(workers)) {
    .sv_positive_integer(workers, "workers")
    if (workers > 1L) {
      warning("workers is retained for compatibility; SVA preprocessing is serial. ",
              "Use blas_threads for a threaded BLAS.")
    }
  }
  1L
}


#### CpG parsing + user set prep

# Accept arbitrary unique feature IDs; recover coordinates when present.
.sv_parse_cpg_ids <- function(cpg_ids) {
  parts <- strsplit(cpg_ids, ".", fixed = TRUE)
  parsed <- lapply(parts, function(p) {
    if (length(p) %in% c(2L, 3L)) {
      x <- suppressWarnings(as.integer(p[2L]))
      y <- if (length(p) == 3L) suppressWarnings(as.integer(p[3L])) else x
      if (!is.na(x) && !is.na(y)) return(list(chr=p[1L], start=x, end=y))
    }
    list(chr=NA_character_, start=NA_integer_, end=NA_integer_)
  })
  data.frame(cpg_id=cpg_ids,
    chr=vapply(parsed, `[[`, character(1), "chr"),
    start=vapply(parsed, `[[`, integer(1), "start"),
    end=vapply(parsed, `[[`, integer(1), "end"), stringsAsFactors=FALSE)
}

.sv_prepare_user_set <- function(all_cpg_ids, user_cpg_set = NULL) {
  if (is.null(user_cpg_set) || length(user_cpg_set) == 0L) {
    return(list(
      user_cpg_set = character(0),
      user_present = character(0),
      n_user_supplied = 0L,
      n_user_present = 0L
    ))
  }
  
  user_cpg_set <- unique(user_cpg_set)
  user_present <- intersect(user_cpg_set, all_cpg_ids)
  
  list(
    user_cpg_set = user_cpg_set,
    user_present = user_present,
    n_user_supplied = length(user_cpg_set),
    n_user_present = length(user_present)
  )
}


#### Missingness + HVF computation

.sv_compute_missingness <- function(meth_mat, use_delayed = TRUE, workers = 1L) {
  if (use_delayed &&
      methods::is(meth_mat, "DelayedMatrix") &&
      requireNamespace("DelayedMatrixStats", quietly = TRUE)) {
    na_count <- DelayedMatrixStats::rowCounts(is.na(meth_mat), value = TRUE)
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    na_count <- matrixStats::rowCounts(is.na(as.matrix(meth_mat)), value = TRUE)
  } else {
    na_count <- rowSums(is.na(as.matrix(meth_mat)))
  }
  
  data.frame(
    cpg_id = rownames(meth_mat),
    n_missing = na_count,
    missing_prop = na_count / ncol(meth_mat),
    stringsAsFactors = FALSE
  )
}

.sv_resolve_n_hvf <- function(n_hvf, override_hvf, n_available) {
  if (n_available < 1L) {
    return(0L)
  }
  
  if (is.null(n_hvf)) {
    if (n_available < 5000L) return(n_available)
    if (n_available < 50000L) return(min(20000L, n_available))
    if (n_available < 200000L) return(min(50000L, n_available))
    return(min(100000L, n_available))
  }
  
  n_hvf <- as.integer(n_hvf)
  
  if (n_hvf < 1L) {
    stop("n_hvf must be >= 1")
  }
  
  if (!override_hvf) {
    if (n_hvf > 150000L) {
      warning("n_hvf > 150000 and override_hvf=FALSE; falling back to 150000")
      n_hvf <- 150000L
    }
  } else {
    if (n_hvf >= n_available) {
      stop("override_hvf=TRUE requires n_hvf < number of available candidate CpGs")
    }
  }
  
  min(n_hvf, n_available)
}

.sv_compute_hvf_stat <- function(meth_mat, candidate_ids, stat = "mad", use_delayed = TRUE) {
  if (length(candidate_ids) == 0L) {
    return(data.frame(cpg_id = character(0), hvf_stat = numeric(0), stringsAsFactors = FALSE))
  }
  
  idx <- match(candidate_ids, rownames(meth_mat))
  x <- meth_mat[idx, , drop = FALSE]
  
  if (use_delayed &&
      methods::is(x, "DelayedMatrix") &&
      requireNamespace("DelayedMatrixStats", quietly = TRUE)) {
    
    if (stat == "mad") {
      vals <- DelayedMatrixStats::rowMads(x, na.rm = TRUE)
    } else {
      vals <- DelayedMatrixStats::rowVars(x, na.rm = TRUE)
    }
    
  } else {
    x_mat <- as.matrix(x)
    # Fast C-level row stats instead of apply() over ~1M rows.
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      if (stat == "mad") {
        vals <- matrixStats::rowMads(x_mat, na.rm = TRUE)
      } else {
        vals <- matrixStats::rowVars(x_mat, na.rm = TRUE)
      }
    } else {
      warning("matrixStats not installed; falling back to slow apply(). install.packages('matrixStats') for a large speedup.")
      if (stat == "mad") {
        vals <- apply(x_mat, 1, stats::mad, na.rm = TRUE)
      } else {
        vals <- apply(x_mat, 1, stats::var, na.rm = TRUE)
      }
    }
  }
  
  vals[is.na(vals)] <- -Inf
  
  data.frame(
    cpg_id = candidate_ids,
    hvf_stat = vals,
    stringsAsFactors = FALSE
  )
}


#### Feature selection

.sv_select_features <- function(
    meth_mat,
    cpg_info,
    user_set_info,
    missing_stats,
    feature_strategy,
    n_hvf,
    override_hvf,
    hvf_stat,
    missingness_max,
    use_delayed,
    workers
) {
  keep_missing <- missing_stats$missing_prop <= missingness_max
  eligible_ids <- missing_stats$cpg_id[keep_missing]
  
  user_ids <- intersect(user_set_info$user_present, eligible_ids)
  n_hvf_used <- 0L
  hvf_tbl <- data.frame(cpg_id = character(0), hvf_stat = numeric(0), stringsAsFactors = FALSE)
  hvf_ids <- character(0)
  
  if (feature_strategy == "user_set_only") {
    selected_ids <- user_ids
    
  } else {
    if (feature_strategy == "hvf_only") {
      hvf_candidates <- eligible_ids
      user_ids <- character(0)
    } else {
      hvf_candidates <- setdiff(eligible_ids, user_ids)
    }
    
    n_hvf_used <- .sv_resolve_n_hvf(
      n_hvf = n_hvf,
      override_hvf = override_hvf,
      n_available = length(hvf_candidates)
    )
    
    if (n_hvf_used > 0L) {
      .svlog(sprintf("Computing %s over %d candidate CpGs", hvf_stat, length(hvf_candidates)))
      hvf_tbl <- .sv_compute_hvf_stat(
        meth_mat = meth_mat,
        candidate_ids = hvf_candidates,
        stat = hvf_stat,
        use_delayed = use_delayed
      )
      
      hvf_tbl <- hvf_tbl[order(hvf_tbl$hvf_stat, decreasing = TRUE), , drop = FALSE]
      hvf_ids <- head(hvf_tbl$cpg_id, n_hvf_used)
    }
    
    selected_ids <- c(user_ids, hvf_ids)
  }
  
  selected_ids <- unique(selected_ids)
  
  feature_table <- data.frame(cpg_id = selected_ids, stringsAsFactors = FALSE)
  feature_table <- merge(feature_table, cpg_info, by = "cpg_id", all.x = TRUE, sort = FALSE)
  feature_table <- merge(feature_table, missing_stats, by = "cpg_id", all.x = TRUE, sort = FALSE)
  
  feature_table$in_user_set <- feature_table$cpg_id %in% user_ids
  feature_table$selected_as_hvf <- feature_table$cpg_id %in% hvf_ids
  
  hvf_map <- hvf_tbl$hvf_stat
  names(hvf_map) <- hvf_tbl$cpg_id
  feature_table$hvf_stat <- unname(hvf_map[feature_table$cpg_id])
  
  feature_table <- feature_table[match(selected_ids, feature_table$cpg_id), , drop = FALSE]
  
  summary <- list(
    n_selected_total = length(selected_ids),
    n_user_set_used = length(user_ids),
    n_hvf_used = length(hvf_ids),
    n_eligible_after_missingness = length(eligible_ids),
    n_user_supplied = user_set_info$n_user_supplied,
    n_user_present = user_set_info$n_user_present
  )
  
  list(
    selected_cpg_ids = selected_ids,
    feature_table = feature_table,
    summary = summary
  )
}


#### Matrix extraction + imputation stub

.sv_extract_selected_matrix <- function(meth_mat, selected_cpg_ids) {
  idx <- match(selected_cpg_ids, rownames(meth_mat))
  x <- meth_mat[idx, , drop = FALSE]
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

.sv_impute_selected_matrix <- function(x, cpg_info, method = "knn_chr", knn_k = 10, workers = 1L) {
  stop(
    paste0(
      "Imputation interface is enabled, but the real chromosome-wise KNN ",
      "imputation engine has not yet been implemented. Set impute=FALSE for now."
    )
  )
}


#### SVA prep + execution

.sv_prepare_matrix_for_sva <- function(x, center_features = TRUE, scale_features = FALSE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  if (center_features || scale_features) {
    x <- t(scale(t(x), center = center_features, scale = scale_features))
  }
  
  keep_finite <- rowSums(!is.finite(x)) == 0L
  if (!all(keep_finite)) {
    x <- x[keep_finite, , drop = FALSE]
  }
  
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    row_var <- matrixStats::rowVars(x)
  } else {
    row_var <- apply(x, 1, stats::var)
  }
  keep_var <- is.finite(row_var) & row_var > 0
  x <- x[keep_var, , drop = FALSE]
  
  if (nrow(x) < 2L) {
    stop("Too few variable CpGs remain after centering/scaling.")
  }
  
  x
}

.sv_choose_num_sv <- function(x, mod, mod0, method = "be", manual_n = NULL) {
  if (method == "manual") {
    return(as.integer(manual_n))
  }
  
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required")
  }
  
  n_sv <- sva::num.sv(dat = x, mod = mod, method = method)
  if (length(n_sv) != 1L || !is.finite(n_sv)) stop("num.sv returned an invalid SV count")
  n_sv <- as.integer(n_sv)

  n_max <- max(0L, min(nrow(x), ncol(x) - qr(mod)$rank - 1L))
  n_sv <- min(n_sv, n_max)
  
  if (ncol(x) < 20L) {
    .svlog("Small sample size detected; n_sv estimate may be unstable")
  }
  
  max(0L, n_sv)
}

.sv_run_sva <- function(x, mod, mod0, n_sv, method = "irw") {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required")
  }
  
  if (method == "irw") {
    fit <- sva::sva(dat = x, mod = mod, mod0 = mod0, n.sv = n_sv)
  } else {
    stop("Use method = 'irw' for methylation M values")
  }
  
  fit
}

.sv_format_sv_output <- function(sva_fit, sample_names) {
  sv <- sva_fit$sv
  sv <- as.data.frame(sv, stringsAsFactors = FALSE)
  colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
  rownames(sv) <- sample_names
  sv$sample_id <- sample_names
  sv <- sv[, c("sample_id", paste0("SV", seq_len(ncol(sva_fit$sv)))), drop = FALSE]
  sv
}


#### Diagnostics

# Headless-safe PNG: compute nodes have no X11 display, so the default
# Xlib png() device fails ("unable to start device PNG"). Use the Cairo
# backend, which renders without X11. Returns TRUE if the device opened.
# On failure it warns and skips the plot rather than killing the job, so
# the accompanying TSV outputs still get written.
.sv_open_png <- function(file, width = 1400, height = 1000, res = 150) {
  ok <- tryCatch({
    if (capabilities("cairo")) {
      grDevices::png(file, width = width, height = height, res = res, type = "cairo")
    } else if (requireNamespace("Cairo", quietly = TRUE)) {
      Cairo::CairoPNG(file, width = width, height = height, res = res)
    } else {
      grDevices::png(file, width = width, height = height, res = res)
    }
    TRUE
  }, error = function(e) {
    warning(sprintf("Could not open PNG device for %s (%s); skipping this plot.",
                    basename(file), conditionMessage(e)))
    FALSE
  })
  ok
}

.sv_generate_diagnostics <- function(
    x,
    sv_df,
    feature_plan,
    cpg_info,
    missing_stats,
    diagnostics,
    diagnostics_dir = NULL
) {
  # Prefer Cairo bitmaps everywhere in this function (belt-and-suspenders
  # alongside type="cairo" in .sv_open_png).
  old_bitmap <- getOption("bitmapType")
  if (isTRUE(capabilities("cairo"))) options(bitmapType = "cairo")
  on.exit({ if (!is.null(old_bitmap)) options(bitmapType = old_bitmap) }, add = TRUE)
  
  out <- list()
  
  if (!is.null(diagnostics_dir)) {
    dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if ("feature_summary" %in% diagnostics) {
    out$feature_summary <- feature_plan$summary
    
    if (!is.null(diagnostics_dir)) {
      utils::write.table(
        as.data.frame(feature_plan$summary, stringsAsFactors = FALSE),
        file = file.path(diagnostics_dir, "feature_summary.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      utils::write.table(
        feature_plan$feature_table,
        file = file.path(diagnostics_dir, "selected_feature_table.tsv.gz"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }
  
  if ("missingness" %in% diagnostics) {
    miss_summary <- summary(missing_stats$missing_prop)
    out$missingness <- miss_summary
    
    if (!is.null(diagnostics_dir)) {
      if (.sv_open_png(file.path(diagnostics_dir, "missingness_hist.png"))) {
        hist(missing_stats$missing_prop, breaks = 50, main = "CpG Missingness", xlab = "Missingness proportion")
        grDevices::dev.off()
      }
      
      utils::write.table(
        data.frame(stat = names(miss_summary), value = as.numeric(miss_summary), stringsAsFactors = FALSE),
        file = file.path(diagnostics_dir, "missingness_summary.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }
  
  if ("sv_cor" %in% diagnostics && ncol(sv_df) > 1L) {
    sv_mat <- as.matrix(sv_df[, setdiff(colnames(sv_df), "sample_id"), drop = FALSE])
    cor_mat <- stats::cor(sv_mat, use = "pairwise.complete.obs")
    out$sv_cor <- cor_mat
    
    if (!is.null(diagnostics_dir)) {
      utils::write.table(
        cor_mat,
        file = file.path(diagnostics_dir, "sv_correlation_matrix.tsv"),
        sep = "\t", quote = FALSE, col.names = NA
      )
    }
  }
  
  if ("scree" %in% diagnostics) {
    scree <- .sv_scree_prop_var(x, n_pc = 20L)
    out$scree <- scree
    
    if (!is.null(diagnostics_dir)) {
      if (.sv_open_png(file.path(diagnostics_dir, "selected_matrix_scree.png"))) {
        plot(
          seq_along(scree), scree, type = "b",
          xlab = "PC",
          ylab = "Proportion variance explained",
          main = "Scree plot of selected matrix"
        )
        grDevices::dev.off()
      }
      
      utils::write.table(
        data.frame(PC = seq_along(scree), prop_var = scree),
        file = file.path(diagnostics_dir, "selected_matrix_scree.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }
  
  if ("sv_pairs" %in% diagnostics && ncol(sv_df) > 2L && !is.null(diagnostics_dir)) {
    sv_cols <- setdiff(colnames(sv_df), "sample_id")
    pairs_to_plot <- utils::combn(sv_cols[seq_len(min(3L, length(sv_cols)))], 2, simplify = FALSE)
    
    for (i in seq_along(pairs_to_plot)) {
      pair <- pairs_to_plot[[i]]
      if (.sv_open_png(
        file.path(diagnostics_dir, sprintf("scatter_%s_vs_%s.png", pair[1], pair[2])),
        width = 1200, height = 1000)) {
        plot(
          sv_df[[pair[1]]], sv_df[[pair[2]]],
          xlab = pair[1], ylab = pair[2],
          main = sprintf("%s vs %s", pair[1], pair[2]),
          pch = 16
        )
        grDevices::dev.off()
      }
    }
  }
  
  if ("runtime" %in% diagnostics && !is.null(diagnostics_dir)) {
    utils::write.table(
      as.data.frame(feature_plan$summary, stringsAsFactors = FALSE),
      file = file.path(diagnostics_dir, "feature_counts_runtime_context.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  
  out
}

# Truncated SVD scree (irlba) instead of a full prcomp of a 200k x N matrix.
.sv_scree_prop_var <- function(x, n_pc = 20L) {
  # Total variance = sum of column variances of the samples-by-features matrix.
  tx <- t(x)  # samples x features
  total_var <- sum(apply(tx, 2, stats::var))
  
  n_pc <- min(n_pc, nrow(tx) - 1L, ncol(tx) - 1L)
  if (n_pc < 1L) n_pc <- 1L
  
  if (requireNamespace("irlba", quietly = TRUE)) {
    pc <- irlba::prcomp_irlba(tx, n = n_pc, center = TRUE, scale. = FALSE)
    prop <- (pc$sdev^2) / total_var
  } else {
    warning("irlba not installed; falling back to full prcomp for scree (slow). install.packages('irlba').")
    pc <- stats::prcomp(tx, center = TRUE, scale. = FALSE)
    prop <- (pc$sdev^2) / sum(pc$sdev^2)
  }
  prop
}


#### Runtime summary

.sv_collect_runtime_summary <- function(start_time, input_dim, x_sel, x_sva, feature_plan, n_sv, workers) {
  end_time <- Sys.time()
  
  list(
    start_time = as.character(start_time),
    end_time = as.character(end_time),
    elapsed_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
    input_dim = input_dim,
    selected_dim = c(nrow(x_sel), ncol(x_sel)),
    sva_dim = c(nrow(x_sva), ncol(x_sva)),
    n_selected_total = feature_plan$summary$n_selected_total,
    n_user_set_used = feature_plan$summary$n_user_set_used,
    n_hvf_used = feature_plan$summary$n_hvf_used,
    n_sv = n_sv,
    workers = workers,
    object_size_selected = as.character(utils::object.size(x_sel)),
    object_size_sva = as.character(utils::object.size(x_sva))
  )
}
