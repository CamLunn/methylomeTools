# =========================================================
# DSS multifactor, block-parallel:
#   - live parallel backend (fork / PSOCK)
#   - slurm array script writer backend
#
# Contiguous within-chromosome blocks balance work when smoothing is disabled.
# Smoothing uses whole chromosomes to avoid artificial block boundaries.
# Input counts must already be filtered according to the analysis plan.
# =========================================================

run_dss_parallel <- function(
    bs,
    meta,
    formula,
    outcome_col,
    outcome_type = c("auto", "continuous", "binary", "categorical"),
    ref = NULL,
    contrast_structure = NULL,
    multi_group_default = c("term", "one_vs_rest", "none"),
    
    # manual override: choose exactly one if using manual mode
    coef = NULL,
    term = NULL,
    Contrast = NULL,
    
    sample_id_col = NULL,
    
    # DSS options
    smoothing = FALSE,
    smoothing.span = 500,
    
    # execution backend
    backend = c("live", "slurm"),
    n_workers = .available_workers(),
    
    parallel_backend = c("psock", "fork"),

    # ---- partitioning / load balancing ----
    partition_mode = c("block", "chromosome"),
    target_block_size = 150000L,   # approx loci per block in block mode
    pin_blas_threads_n = 1L,       # threads per worker for BLAS/OpenMP
    
    # chromosome selection
    chr_subset = NULL,
    canonical_chr_only = TRUE,
    
    # outputs
    outdir = file.path(getwd(), "dss_multifactor_by_chr"),
    result_prefix = "dss_mf",
    save_rds = TRUE,
    save_fit = FALSE,
    save_compression = c("gzip", "none", "xz"),
    save_test_only = TRUE,
    save_log_memory = TRUE,
    auto_merge = TRUE,
    resume = FALSE,
    
    # slurm options
    slurm_partition = NULL,
    slurm_time = "24:00:00",
    slurm_mem = "16G",             # per-task; block subsets are small now
    slurm_cpus_per_task = 1,
    slurm_account = NULL,
    slurm_module_load = character(),
    rscript_path = "Rscript",
    slurm_subset_compress = "gzip",  # compression for per-block subset objects
    
    # logging
    verbose = TRUE,
    gc_each_chr = TRUE
) {
  outcome_type <- match.arg(outcome_type)
  multi_group_default <- match.arg(multi_group_default)
  backend <- match.arg(backend)
  parallel_backend <- match.arg(parallel_backend)
  n_workers <- as.integer(n_workers)
  if (length(n_workers) != 1L || is.na(n_workers) || n_workers < 1L)
    stop("n_workers must be a positive integer")
  partition_mode <- match.arg(partition_mode)
  
  save_compression <- match.arg(save_compression)
  compress_opt <- switch(
    save_compression,
    gzip = "gzip",
    xz   = "xz",
    none = FALSE
  )
  
  msg <- function(...) {
    if (isTRUE(verbose)) message(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  }
  
  # ----------------------------
  # BLAS / OpenMP thread pinning
  # Best effort in the parent; each worker also pins itself.
  # ----------------------------
  pin_blas_threads <- function(n = 1L) {
    n <- as.character(as.integer(n))
    Sys.setenv(
      OMP_NUM_THREADS = n,
      OPENBLAS_NUM_THREADS = n,
      MKL_NUM_THREADS = n,
      VECLIB_MAXIMUM_THREADS = n,
      NUMEXPR_NUM_THREADS = n
    )
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(as.integer(n)), silent = TRUE)
      try(RhpcBLASctl::omp_set_num_threads(as.integer(n)),  silent = TRUE)
    }
    invisible(NULL)
  }
  pin_blas_threads(pin_blas_threads_n)
  
  # block mode + smoothing is unsafe (boundary edge effects). Fall back.
  if (isTRUE(smoothing) && partition_mode == "block") {
    warning("smoothing = TRUE is incompatible with block partitioning ",
            "(boundary edge effects). Falling back to partition_mode = 'chromosome'.")
    partition_mode <- "chromosome"
  }
  
  # ----------------------------
  # package checks
  # ----------------------------
  req_pkgs <- c("DSS", "bsseq", "GenomicRanges", "SummarizedExperiment")
  for (pkg in req_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Missing package: ", pkg)
    }
  }
  
  # ----------------------------
  # directory layout
  # ----------------------------
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  run_config <- list(samples = colnames(bs), loci = SummarizedExperiment::rowRanges(bs),
    formula = paste(deparse(formula), collapse=" "), meta = meta, outcome_col=outcome_col,
    outcome_type=outcome_type, ref=ref, coef=coef, term=term, Contrast=Contrast,
    contrast_structure=contrast_structure, multi_group_default=multi_group_default,
    smoothing=smoothing, smoothing.span=smoothing.span, partition_mode=partition_mode,
    target_block_size=target_block_size, chr_subset=chr_subset, canonical_chr_only=canonical_chr_only)
  config_file <- file.path(outdir, paste0(result_prefix, "_run_config.rds"))
  if (!isTRUE(resume) && dir.exists(file.path(outdir, "results_by_chr")) &&
      length(list.files(file.path(outdir, "results_by_chr"), pattern="\\.rds$")))
    stop("Results already exist. Choose a new outdir or explicitly set resume=TRUE for unchanged inputs.")
  if (file.exists(config_file) && !identical(readRDS(config_file), run_config))
    stop("Output directory contains a different run configuration. Choose a new outdir or result_prefix.")
  if (!file.exists(config_file)) {
    if (dir.exists(file.path(outdir, "results_by_chr")) &&
        length(list.files(file.path(outdir, "results_by_chr"), pattern="\\.rds$")))
      stop("Existing results have no run configuration. Choose a new output directory.")
    saveRDS(run_config, config_file)
  }
  logs_dir    <- file.path(outdir, "logs")
  results_dir <- file.path(outdir, "results_by_chr")
  scripts_dir <- file.path(outdir, "scripts")
  meta_dir    <- file.path(outdir, "meta")
  units_dir   <- file.path(outdir, "unit_inputs")  # per-block subset objects (slurm)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(units_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------
  # helpers
  # ----------------------------
  get_bs_samples <- function(bs_obj) {
    out <- tryCatch(bsseq::sampleNames(bs_obj), error = function(e) NULL)
    if (is.null(out) || length(out) == 0) out <- colnames(bs_obj)
    out
  }
  
  get_chr_levels <- function(bs_obj) {
    rr <- SummarizedExperiment::rowRanges(bs_obj)
    chrs <- as.character(GenomicRanges::seqnames(rr))
    unique(chrs)
  }
  
  is_canonical_chr <- function(chr) {
    grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", chr)
  }
  
  normalize_char_levels <- function(x) {
    trimws(as.character(x))
  }
  
  coerce_outcome_type_auto <- function(x) {
    x_non_na <- x[!is.na(x)]
    if (is.numeric(x_non_na) || is.integer(x_non_na)) return("continuous")
    if (is.logical(x_non_na)) return("binary")
    if (is.factor(x_non_na) || is.character(x_non_na)) {
      u <- unique(as.character(x_non_na))
      if (length(u) == 2) return("binary")
      return("categorical")
    }
    "categorical"
  }
  
  validate_contrast_structure <- function(cs, levels_outcome) {
    if (!is.list(cs) || is.null(names(cs)) || any(names(cs) == "")) {
      stop("contrast_structure must be a named list.")
    }
    for (nm in names(cs)) {
      one <- cs[[nm]]
      if (!is.list(one) || !all(c("case", "control") %in% names(one))) {
        stop("Each contrast entry must be a list with case and control: ", nm)
      }
      case <- normalize_char_levels(one$case)
      control <- normalize_char_levels(one$control)
      if (length(case) < 1 || length(control) < 1) {
        stop("Contrast '", nm, "' has empty case/control.")
      }
      if (any(!case %in% levels_outcome)) {
        stop("Unknown case levels in contrast '", nm, "': ",
             paste(setdiff(case, levels_outcome), collapse = ", "))
      }
      if (any(!control %in% levels_outcome)) {
        stop("Unknown control levels in contrast '", nm, "': ",
             paste(setdiff(control, levels_outcome), collapse = ", "))
      }
      if (length(intersect(case, control)) > 0) {
        stop("Overlap between case/control in contrast '", nm, "'.")
      }
    }
    invisible(TRUE)
  }
  
  build_binary_from_levels <- function(x, case_levels, control_levels) {
    grp <- ifelse(as.character(x) %in% case_levels, "case",
                  ifelse(as.character(x) %in% control_levels, "control", NA_character_))
    factor(grp, levels = c("control", "case"))
  }
  
  get_memory_snapshot_mb <- function() {
    rss_mb <- NA_real_
    vmhwm_mb <- NA_real_
    
    status_file <- "/proc/self/status"
    if (file.exists(status_file)) {
      txt <- tryCatch(readLines(status_file, warn = FALSE), error = function(e) character())
      if (length(txt)) {
        rss_line <- grep("^VmRSS:", txt, value = TRUE)
        hwm_line <- grep("^VmHWM:", txt, value = TRUE)
        parse_kb <- function(z) {
          if (!length(z)) return(NA_real_)
          as.numeric(gsub("[^0-9]", "", z)) / 1024
        }
        rss_mb <- parse_kb(rss_line)
        vmhwm_mb <- parse_kb(hwm_line)
      }
    }
    
    gc_mat <- gc()
    gc_used_mb <- sum(gc_mat[, "used"] * c(8, 8)) / 1024^2
    
    data.frame(
      rss_mb = rss_mb,
      vmhwm_mb = vmhwm_mb,
      gc_used_mb = gc_used_mb,
      stringsAsFactors = FALSE
    )
  }
  
  write_log_line <- function(log_file, ...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
    cat(line, "\n", file = log_file, append = TRUE)
  }
  
  safe_save_rds <- function(obj, path, compress = "gzip") {
    saveRDS(obj, path, compress = compress)
    invisible(path)
  }
  
  # ----------------------------
  # metadata alignment + checks
  # ----------------------------
  samp <- get_bs_samples(bs)
  if (is.null(samp) || length(samp) == 0) {
    stop("Could not determine sample names from BSseq object.")
  }
  
  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(meta)) {
      stop("sample_id_col not found in meta: ", sample_id_col)
    }
    rownames(meta) <- as.character(meta[[sample_id_col]])
  }
  
  if (is.null(rownames(meta))) {
    stop("meta must have rownames matching BSseq sample names, or provide sample_id_col.")
  }
  
  if (anyDuplicated(rownames(meta))) {
    dup <- unique(rownames(meta)[duplicated(rownames(meta))])
    stop("Duplicated sample IDs in meta: ", paste(dup, collapse = ", "))
  }
  
  missing_samples <- setdiff(samp, rownames(meta))
  if (length(missing_samples)) {
    stop("meta missing BSseq samples: ", paste(missing_samples, collapse = ", "))
  }
  
  meta2 <- droplevels(meta[samp, , drop = FALSE])
  
  if (is.character(formula)) formula <- stats::as.formula(formula)
  if (!inherits(formula, "formula")) {
    stop("formula must be a formula or a character convertible to formula.")
  }
  if (length(formula) != 2) {
    stop("formula must be one-sided, e.g. ~ outcome + batch + age")
  }
  
  if (!outcome_col %in% colnames(meta2)) {
    stop("outcome_col not found in aligned meta: ", outcome_col)
  }
  
  vars_needed <- unique(all.vars(formula))
  vars_needed <- vars_needed[vars_needed %in% colnames(meta2)]
  
  na_counts <- vapply(vars_needed, function(v) sum(is.na(meta2[[v]])), numeric(1))
  bad_na <- names(na_counts)[na_counts > 0]
  if (length(bad_na)) {
    stop(
      "Missing values found in metadata model columns: ",
      paste(bad_na, collapse = ", "),
      ".\nThis wrapper requires complete metadata columns for the design matrix."
    )
  }
  
  X_preview <- model.matrix(formula, data = meta2)
  if (qr(X_preview)$rank < ncol(X_preview)) {
    stop("Design matrix is rank deficient before DSS fitting. Check formula/covariates.")
  }
  
  # determine manual mode
  n_manual <- sum(!is.null(coef), !is.null(term), !is.null(Contrast))
  manual_mode <- n_manual > 0
  if (manual_mode && n_manual != 1) {
    stop("Provide exactly one of coef, term, or Contrast in manual mode.")
  }
  
  outcome_type_final <- if (outcome_type == "auto") {
    coerce_outcome_type_auto(meta2[[outcome_col]])
  } else {
    outcome_type
  }
  
  # ----------------------------
  # chromosome inventory
  # ----------------------------
  chr_all <- get_chr_levels(bs)
  if (canonical_chr_only) {
    chr_all <- chr_all[is_canonical_chr(chr_all)]
  }
  if (!is.null(chr_subset)) {
    chr_all <- intersect(chr_all, chr_subset)
  }
  chr_all <- unique(chr_all)
  
  if (!length(chr_all)) {
    stop("No chromosomes selected.")
  }
  
  # ----------------------------
  # work-unit construction (blocks or whole chromosomes)
  # A unit = list(label, chr, idx) where idx are global row indices into bs.
  # ----------------------------
  build_work_units <- function(bs_obj, chr_keep, mode, block_size) {
    rr <- SummarizedExperiment::rowRanges(bs_obj)
    chr_vec <- as.character(GenomicRanges::seqnames(rr))
    units <- list()
    for (this_chr in chr_keep) {
      ci <- which(chr_vec == this_chr)            # position-sorted global indices
      if (!length(ci)) next
      if (mode == "chromosome") {
        units[[length(units) + 1L]] <- list(label = this_chr, chr = this_chr, idx = ci)
      } else {
        grp <- ceiling(seq_along(ci) / block_size)  # contiguous fixed-size blocks
        splits <- split(ci, grp)
        nb <- length(splits)
        for (b in seq_len(nb)) {
          lab <- if (nb == 1L) this_chr else sprintf("%s_b%02d", this_chr, b)
          units[[length(units) + 1L]] <- list(label = lab, chr = this_chr, idx = splits[[b]])
        }
      }
    }
    if (!length(units)) stop("No work units created.")
    # longest-processing-time-first: hand out biggest blocks earliest
    ord <- order(vapply(units, function(u) length(u$idx), integer(1)), decreasing = TRUE)
    units[ord]
  }
  
  work_units <- build_work_units(bs, chr_all, partition_mode, as.integer(target_block_size))
  
  unit_manifest <- data.frame(
    task_id = seq_along(work_units),
    label   = vapply(work_units, function(u) u$label, character(1)),
    chr     = vapply(work_units, function(u) u$chr,   character(1)),
    n_loci  = vapply(work_units, function(u) length(u$idx), integer(1)),
    stringsAsFactors = FALSE
  )
  
  msg("Partition mode: ", partition_mode,
      " | chromosomes: ", length(chr_all),
      " | work units: ", nrow(unit_manifest),
      " | loci/unit (min/median/max): ",
      paste(stats::quantile(unit_manifest$n_loci, c(0, .5, 1)), collapse = "/"))
  
  manifest_file <- file.path(meta_dir, paste0(result_prefix, "_work_manifest.tsv"))
  write.table(unit_manifest, manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  meta_file <- file.path(meta_dir, paste0(result_prefix, "_meta.rds"))
  bs_file   <- file.path(meta_dir, paste0(result_prefix, "_bs.rds"))
  
  # ----------------------------
  # inner test planner
  # ----------------------------
  build_test_plan <- function(meta_df, formula, outcome_col, outcome_type_final, ref,
                              contrast_structure, multi_group_default,
                              coef, term, Contrast) {
    if (manual_mode) {
      return(list(
        mode = "manual",
        test_plan = list(
          list(
            name = "manual",
            type = if (!is.null(coef)) "coef" else if (!is.null(term)) "term" else "contrast",
            coef = coef,
            term = term,
            Contrast = Contrast
          )
        ),
        meta = meta_df,
        formula = formula,
        outcome_type_final = outcome_type_final
      ))
    }
    
    term_labels <- attr(stats::terms(formula), "term.labels")
    cov_terms   <- setdiff(term_labels, outcome_col)
    
    # ---- CONTINUOUS OUTCOME ----
    if (outcome_type_final == "continuous") {
      formula_fit <- stats::reformulate(c(outcome_col, cov_terms))
      X_preview <- model.matrix(formula_fit, data = meta_df)
      
      coef_name <- outcome_col
      if (!outcome_col %in% colnames(X_preview)) {
        hits <- grep(
          paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", outcome_col)),
          colnames(X_preview)
        )
        if (length(hits) == 1) {
          coef_name <- colnames(X_preview)[hits]
        } else if (length(hits) > 1) {
          stop("Ambiguous design columns for outcome '", outcome_col,
               "': ", paste(colnames(X_preview)[hits], collapse = ", "),
               ". Please specify `coef` manually.")
        } else {
          stop("Could not find a design matrix column corresponding to outcome '",
               outcome_col, "'. Design columns: ",
               paste(colnames(X_preview), collapse = ", "))
        }
      }
      
      return(list(
        mode = "auto",
        test_plan = list(list(name = "main", type = "coef", coef = coef_name)),
        meta = meta_df,
        formula = formula_fit,
        outcome_type_final = outcome_type_final
      ))
    }
    
    # ---- BINARY OUTCOME ----
    if (outcome_type_final == "binary") {
      x <- meta_df[[outcome_col]]
      if (is.logical(x)) x <- ifelse(x, "TRUE", "FALSE")
      x <- factor(normalize_char_levels(x))
      lv <- levels(x)
      if (length(lv) != 2) stop("Binary outcome requires exactly 2 levels after coercion.")
      if (is.null(ref) || !ref %in% lv) ref <- lv[1]
      x <- stats::relevel(x, ref = ref)
      meta_df[[outcome_col]] <- x
      
      formula_fit <- stats::reformulate(c(outcome_col, cov_terms))
      X_preview <- model.matrix(formula_fit, data = meta_df)
      
      non_ref_level <- setdiff(levels(x), ref)
      if (length(non_ref_level) != 1) {
        stop("Internal error: binary outcome should have exactly one non-reference level. ",
             "Got: ", paste(non_ref_level, collapse = ", "))
      }
      lvl <- non_ref_level
      cand1 <- paste0(outcome_col, lvl)
      cand2 <- paste0(outcome_col, make.names(lvl))
      coef_idx <- which(colnames(X_preview) %in% c(cand1, cand2))
      
      if (length(coef_idx) == 1) {
        coef_name <- colnames(X_preview)[coef_idx]
      } else {
        hits <- grep(
          paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", outcome_col)),
          colnames(X_preview)
        )
        hits <- setdiff(hits, match("(Intercept)", colnames(X_preview), nomatch = 0))
        if (length(hits) == 1) {
          coef_name <- colnames(X_preview)[hits]
        } else if (length(hits) > 1) {
          stop("Ambiguous design columns for binary outcome '", outcome_col,
               "': ", paste(colnames(X_preview)[hits], collapse = ", "),
               ". Please specify `coef` manually.")
        } else {
          stop("Could not find a design matrix column corresponding to binary outcome '",
               outcome_col, "'. Design columns: ",
               paste(colnames(X_preview), collapse = ", "))
        }
      }
      
      return(list(
        mode = "auto",
        test_plan = list(list(name = "main", type = "coef", coef = coef_name)),
        meta = meta_df,
        formula = formula_fit,
        outcome_type_final = outcome_type_final
      ))
    }
    
    # ---- CATEGORICAL OUTCOME ----
    if (outcome_type_final == "categorical") {
      x <- factor(normalize_char_levels(meta_df[[outcome_col]]))
      lv <- levels(x)
      if (length(lv) < 3) stop("Categorical outcome requires 3+ levels.")
      meta_df[[outcome_col]] <- x
      
      formula_fit <- stats::reformulate(c(outcome_col, cov_terms))
      plan <- list()
      
      if (multi_group_default == "term") {
        plan[[length(plan) + 1L]] <- list(
          name = paste0(outcome_col, "_omnibus"), type = "term", term = outcome_col
        )
      }
      
      if (multi_group_default == "one_vs_rest") {
        if (is.null(contrast_structure)) {
          contrast_structure <- setNames(
            lapply(lv, function(L) list(case = L, control = setdiff(lv, L))),
            paste0(lv, "_vs_rest")
          )
        }
      }
      
      if (!is.null(contrast_structure)) {
        validate_contrast_structure(contrast_structure, lv)
        for (nm in names(contrast_structure)) {
          plan[[length(plan) + 1L]] <- list(
            name = nm, type = "level_contrast",
            case = normalize_char_levels(contrast_structure[[nm]]$case),
            control = normalize_char_levels(contrast_structure[[nm]]$control),
            outcome_col = outcome_col
          )
        }
      }
      
      if (!length(plan)) stop("No categorical test plan created. Use multi_group_default or contrast_structure.")
      
      return(list(
        mode = "auto",
        test_plan = plan,
        meta = meta_df,
        formula = formula_fit,
        outcome_type_final = outcome_type_final
      ))
    }
    
    stop("Unhandled outcome_type_final: ", outcome_type_final)
  }
  
  # ----------------------------
  # build Contrast matrix for factor level comparisons
  # ----------------------------
  build_factor_level_contrast <- function(fit, outcome_col, case, control, meta_df) {
    xcols <- colnames(fit$X)
    levs <- levels(meta_df[[outcome_col]])
    ref_level <- levs[1]
    
    coef_for_level <- function(level_name) {
      if (identical(level_name, ref_level)) return(NULL)
      exact1 <- paste0(outcome_col, level_name)
      exact2 <- paste0(outcome_col, make.names(level_name))
      hits <- which(xcols %in% c(exact1, exact2))
      if (length(hits) == 1) return(hits)
      hits <- grep(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", outcome_col)), xcols)
      hit_names <- xcols[hits]
      hit_idx <- hits[grepl(make.names(level_name), make.names(hit_names), fixed = TRUE) |
                        grepl(level_name, hit_names, fixed = TRUE)]
      if (length(hit_idx) == 1) return(hit_idx)
      stop("Could not uniquely map factor level '", level_name,
           "' to fit$X columns. Available: ", paste(xcols, collapse = ", "))
    }
    
    v <- rep(0, ncol(fit$X))
    w_case <- 1 / length(case)
    w_ctrl <- 1 / length(control)
    for (lv in case) { idx <- coef_for_level(lv); if (!is.null(idx)) v[idx] <- v[idx] + w_case }
    for (lv in control) { idx <- coef_for_level(lv); if (!is.null(idx)) v[idx] <- v[idx] - w_ctrl }
    matrix(v, ncol = 1)
  }
  
  # ----------------------------
  # run one work unit (block or chromosome)
  # ----------------------------
  run_one_unit <- function(unit, bs_obj, meta_df, formula, outcome_col,
                           outcome_type_final, ref, contrast_structure,
                           multi_group_default, coef, term, Contrast,
                           smoothing, smoothing.span,
                           results_dir, logs_dir,
                           result_prefix, save_rds, save_fit,
                           save_log_memory = TRUE,
                           save_compression = "gzip",
                           gc_each_chr = FALSE,
                           pin_blas_threads_n = 1L) {
    
    # each worker pins its own threads (robust for forked / PSOCK procs)
    Sys.setenv(OMP_NUM_THREADS = pin_blas_threads_n,
               OPENBLAS_NUM_THREADS = pin_blas_threads_n,
               MKL_NUM_THREADS = pin_blas_threads_n,
               VECLIB_MAXIMUM_THREADS = pin_blas_threads_n)
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(as.integer(pin_blas_threads_n)), silent = TRUE)
      try(RhpcBLASctl::omp_set_num_threads(as.integer(pin_blas_threads_n)),  silent = TRUE)
    }
    
    label  <- unit$label
    chr    <- unit$chr
    idx    <- unit$idx
    
    log_file <- file.path(logs_dir, paste0(result_prefix, "_", label, ".log"))
    res_file <- file.path(results_dir, paste0(result_prefix, "_", label, ".rds"))
    
    compress_opt <- switch(save_compression, gzip = "gzip", xz = "xz", none = FALSE, "gzip")
    
    t0_total <- proc.time()[["elapsed"]]
    mem0 <- get_memory_snapshot_mb()
    write_log_line(log_file, "Starting unit: ", label, " (chr=", chr, ")")
    write_log_line(log_file, "Initial memory (MB): rss=", round(mem0$rss_mb, 1),
                   ", vmhwm=", round(mem0$vmhwm_mb, 1),
                   ", gc_used=", round(mem0$gc_used_mb, 1))
    
    n_loci_chr <- length(idx)
    if (n_loci_chr == 0) {
      out <- list(chr = chr, label = label, status = "empty", n_loci = 0L, result_file = NA_character_)
      if (save_rds) safe_save_rds(out, res_file, compress = compress_opt)
      write_log_line(log_file, "No loci in unit.")
      return(out)
    }
    
    write_log_line(log_file, "Loci in unit: ", n_loci_chr)
    bs_chr <- bs_obj[idx, ]   # subset by precomputed indices (block or whole chr)
    
    test_plan_obj <- build_test_plan(
      meta_df = meta_df, formula = formula, outcome_col = outcome_col,
      outcome_type_final = outcome_type_final, ref = ref,
      contrast_structure = contrast_structure, multi_group_default = multi_group_default,
      coef = coef, term = term, Contrast = Contrast
    )
    
    meta_use <- test_plan_obj$meta
    formula_use <- test_plan_obj$formula
    
    X_preview <- model.matrix(formula_use, data = meta_use)
    if (qr(X_preview)$rank < ncol(X_preview)) {
      stop("Rank-deficient design matrix on ", label)
    }
    
    write_log_line(log_file, "Formula used: ", paste(format(formula_use), collapse = " "))
    write_log_line(log_file, "Design columns preview: ", paste(colnames(X_preview), collapse = ", "))
    
    t0_fit <- proc.time()[["elapsed"]]
    fit <- DSS::DMLfit.multiFactor(
      BSobj = bs_chr, design = meta_use, formula = formula_use,
      smoothing = smoothing, smoothing.span = smoothing.span
    )
    t_fit <- proc.time()[["elapsed"]] - t0_fit
    
    write_log_line(log_file, "Fit complete. fit$X cols: ", paste(colnames(fit$X), collapse = ", "))
    write_log_line(log_file, "Fit time (sec): ", round(t_fit, 2))
    
    tests_out <- list()
    for (tp in test_plan_obj$test_plan) {
      write_log_line(log_file, "Running test: ", tp$name, " [", tp$type, "]")
      t0_test <- proc.time()[["elapsed"]]
      
      test_res <- switch(
        tp$type,
        coef = DSS::DMLtest.multiFactor(fit, coef = tp$coef),
        term = DSS::DMLtest.multiFactor(fit, term = tp$term),
        contrast = DSS::DMLtest.multiFactor(fit, Contrast = tp$Contrast),
        level_contrast = {
          Cmat <- build_factor_level_contrast(
            fit = fit, outcome_col = tp$outcome_col,
            case = tp$case, control = tp$control, meta_df = meta_use
          )
          DSS::DMLtest.multiFactor(fit, Contrast = Cmat)
        },
        stop("Unknown test type: ", tp$type)
      )
      
      t_test <- proc.time()[["elapsed"]] - t0_test
      write_log_line(log_file, "Test done: ", tp$name,
                     " | loci=", nrow(test_res), " | sec=", round(t_test, 2))
      
      tests_out[[tp$name]] <- list(
        test = test_res,
        summary = list(name = tp$name, type = tp$type, time_sec = t_test, n_loci = nrow(test_res))
      )
    }
    
    mem1 <- get_memory_snapshot_mb()
    t_total_pre_save <- proc.time()[["elapsed"]] - t0_total
    
    compact_tests <- lapply(tests_out, function(x) x$test)
    compact_summaries <- lapply(tests_out, function(x) x$summary)
    
    out_small <- list(
      chr = chr,
      label = label,
      status = "ok",
      n_loci = n_loci_chr,
      tests = compact_tests,
      test_summaries = compact_summaries,
      timing_sec = list(fit = t_fit, pre_save_total = t_total_pre_save)
    )
    
    if (isTRUE(save_fit)) out_small$fit <- fit
    if (isTRUE(save_log_memory)) out_small$memory_mb <- list(start = mem0, end = mem1)
    
    t0_save <- proc.time()[["elapsed"]]
    write_log_line(log_file, "Beginning save: ", res_file, " | compression=", as.character(compress_opt))
    if (save_rds) safe_save_rds(out_small, res_file, compress = compress_opt)
    t_save <- proc.time()[["elapsed"]] - t0_save
    t_total <- proc.time()[["elapsed"]] - t0_total
    
    write_log_line(log_file, "Finished save in sec: ", round(t_save, 2))
    write_log_line(log_file, "Finished unit ", label,
                   " | total sec=", round(t_total, 2),
                   " | final rss=", round(mem1$rss_mb, 1),
                   " | final vmhwm=", round(mem1$vmhwm_mb, 1))
    
    if (isTRUE(gc_each_chr)) {
      rm(bs_chr, fit, tests_out, out_small, compact_tests, compact_summaries)
      invisible(gc())
      write_log_line(log_file, "Ran gc() after unit completion.")
    }
    
    list(
      chr = chr, label = label, status = "ok", n_loci = n_loci_chr,
      result_file = if (save_rds) res_file else NA_character_,
      log_file = log_file, fit_time_sec = t_fit, save_time_sec = t_save,
      total_time_sec = t_total, rss_end_mb = mem1$rss_mb, vmhwm_end_mb = mem1$vmhwm_mb
    )
  }
  
  recalculate_global_fdr <- function(
    merged_df, p_col = "pvals", old_fdr_col = "fdrs",
    new_fdr_col = "fdr_global", overwrite_fdr = FALSE
  ) {
    if (!is.data.frame(merged_df)) stop("merged_df must be a data.frame")
    if (!p_col %in% colnames(merged_df)) stop("Column not found: ", p_col)
    
    p <- merged_df[[p_col]]
    if (!is.numeric(p)) stop("`", p_col, "` must be numeric.")
    
    if (!is.null(old_fdr_col) && old_fdr_col %in% colnames(merged_df)) {
      backup_col <- paste0(old_fdr_col, "_per_chr")
      if (!backup_col %in% colnames(merged_df)) {
        merged_df[[backup_col]] <- merged_df[[old_fdr_col]]
      }
    }
    
    merged_df[[new_fdr_col]] <- stats::p.adjust(p, method = "BH")
    if (isTRUE(overwrite_fdr)) merged_df[[old_fdr_col]] <- merged_df[[new_fdr_col]]
    merged_df
  }
  
  merge_dss_results <- function(
    results_dir, result_prefix, merged_dir = results_dir,
    p_col = "pvals", old_fdr_col = "fdrs", new_fdr_col = "fdr_global",
    overwrite_fdr = FALSE, save_merged_rds = TRUE, save_merged_tsv = TRUE
  ) {
    files <- list.files(results_dir,
                        pattern = paste0("^", result_prefix, "_.*\\.rds$"),
                        full.names = TRUE)
    if (!length(files)) stop("No per-unit result files found.")
    
    res_list <- lapply(files, readRDS)
    res_ok <- Filter(function(x) identical(x$status, "ok"), res_list)
    if (!length(res_ok)) stop("No successful unit result files found.")
    
    shared_tests <- Reduce(intersect, lapply(res_ok, function(x) names(x$tests)))
    if (!length(shared_tests)) stop("No shared test names across unit files.")
    
    merged_tables <- list()
    for (nm in shared_tests) {
      tabs <- lapply(res_ok, function(x) {
        tab <- x$tests[[nm]]
        if (is.null(tab)) return(NULL)
        tab$chr_job <- x$chr
        tab
      })
      tabs <- Filter(Negate(is.null), tabs)
      merged_tab <- do.call(rbind, tabs)
      rownames(merged_tab) <- NULL
      
      merged_tab <- recalculate_global_fdr(
        merged_df = merged_tab, p_col = p_col, old_fdr_col = old_fdr_col,
        new_fdr_col = new_fdr_col, overwrite_fdr = overwrite_fdr
      )
      merged_tables[[nm]] <- merged_tab
      
      if (isTRUE(save_merged_tsv)) {
        out_tsv <- file.path(merged_dir, paste0(result_prefix, "_", nm, "_merged.tsv.gz"))
        con <- gzfile(out_tsv, open = "wt")
        utils::write.table(merged_tab, file = con, sep = "\t", row.names = FALSE, quote = FALSE)
        close(con)
      }
    }
    
    merged_file <- file.path(merged_dir, paste0(result_prefix, "_merged_results.rds"))
    if (isTRUE(save_merged_rds)) saveRDS(merged_tables, merged_file, compress = "gzip")
    
    summary_df <- do.call(rbind, lapply(res_ok, function(x) {
      data.frame(chr = x$chr, label = x$label %||% x$chr, n_loci = x$n_loci, stringsAsFactors = FALSE)
    }))
    summary_file <- file.path(merged_dir, paste0(result_prefix, "_merged_summary.tsv"))
    utils::write.table(summary_df, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    list(merged_tables = merged_tables,
         merged_file = if (isTRUE(save_merged_rds)) merged_file else NA_character_,
         summary_file = summary_file)
  }
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # ----------------------------
  # live backend
  # ----------------------------
  run_live_backend <- function() {
    msg("Live backend: ", n_workers, " workers over ", length(work_units), " work units.")
    msg("Output directory: ", outdir)
    
    run_fun <- function(unit) {
      res_file <- file.path(results_dir, paste0(result_prefix, "_", unit$label, ".rds"))
      if (file.exists(res_file) &&
          identical(tryCatch(readRDS(res_file)$status, error = function(e) NA), "ok")) {
        return(list(chr = unit$chr, label = unit$label, status = "ok",
                    n_loci = NA, result_file = res_file))   # already done, skip
      }
      tryCatch(
        run_one_unit(
          unit = unit, bs_obj = bs, meta_df = meta2, formula = formula,
          outcome_col = outcome_col, outcome_type_final = outcome_type_final,
          ref = ref, contrast_structure = contrast_structure,
          multi_group_default = multi_group_default, coef = coef, term = term,
          Contrast = Contrast, smoothing = smoothing, smoothing.span = smoothing.span,
          results_dir = results_dir, logs_dir = logs_dir, result_prefix = result_prefix,
          save_rds = save_rds, save_fit = save_fit,
          save_log_memory = save_log_memory, save_compression = save_compression,
          gc_each_chr = gc_each_chr,
          pin_blas_threads_n = pin_blas_threads_n
        ),
        error = function(e) {
          log_file <- file.path(logs_dir, paste0(result_prefix, "_", unit$label, ".log"))
          write_log_line(log_file, "ERROR: ", conditionMessage(e))
          list(chr = unit$chr, label = unit$label, status = "error",
               error = conditionMessage(e), log_file = log_file)
        }
      )
    }
    
    if (n_workers == 1L) {
      res <- lapply(work_units, run_fun)
    } else if (parallel_backend == "fork" && .Platform$OS.type == "unix") {
      # dynamic dispatch: largest blocks (front of work_units) handed out first,
      # smaller blocks backfill as cores free -> balanced wall-clock.
      res <- parallel::mclapply(
        work_units, FUN = run_fun,
        mc.cores = n_workers, mc.preschedule = FALSE
      )
    } else {
      cl <- parallel::makeCluster(n_workers)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      # Fresh R sessions need the same libraries and S4 methods as the parent.
      parallel::clusterCall(cl, function(paths) {
        .libPaths(paths)
        for (pkg in c("methylomeTools", "DSS", "bsseq", "GenomicRanges",
                      "SummarizedExperiment")) {
          loadNamespace(pkg)
        }
        NULL
      }, .libPaths())
      parallel::clusterExport(
        cl,
        varlist = c(
          "bs", "meta2", "formula", "outcome_col", "outcome_type_final",
          "ref", "contrast_structure", "multi_group_default", "coef", "term", "Contrast",
          "smoothing", "smoothing.span", "results_dir", "logs_dir", "result_prefix",
          "save_rds", "save_fit", "save_log_memory", "save_compression",
          "gc_each_chr", "pin_blas_threads_n",
          "run_one_unit", "build_test_plan", "build_factor_level_contrast",
          "normalize_char_levels", "validate_contrast_structure",
          "build_binary_from_levels", "get_memory_snapshot_mb",
          "write_log_line", "safe_save_rds"
        ),
        envir = environment()
      )
      parallel::clusterEvalQ(cl, {
        Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      })
      res <- parallel::parLapplyLB(cl, work_units, run_fun)  # load-balanced PSOCK
    }
    
    # Worker failures and resumed units need not have identical summary fields.
    res <- lapply(seq_along(work_units), function(i) {
      item <- res[[i]]
      if (!is.list(item) || is.null(item$status))
        return(list(chr=work_units[[i]]$chr, label=work_units[[i]]$label,
                    status="error", error="Worker returned no valid result"))
      item
    })
    res <- lapply(res, function(item) {
      if (!is.null(item$error)) item$error <- gsub("[\r\n\t]+", " ", item$error)
      item
    })
    res_df <- dplyr::bind_rows(res)
    summary_file <- file.path(outdir, paste0(result_prefix, "_live_summary.tsv"))
    write.table(res_df, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
    msg("Live run complete. Summary: ", summary_file)
    
    failed <- vapply(res, function(item) !identical(item$status, "ok"), logical(1))
    if (any(failed)) {
      details <- vapply(res[failed], function(item) {
        reason <- if (is.null(item$error)) "No successful result returned" else item$error
        log_path <- if (is.null(item$log_file)) "" else paste0(" [", item$log_file, "]")
        paste0(item$label, ": ", reason, log_path)
      }, character(1))
      stop("Incomplete run: ", sum(failed), " of ", length(res),
           " work units failed.\n", paste(details, collapse = "\n"),
           "\nSummary: ", summary_file, call. = FALSE)
    }
    if (isTRUE(save_rds)) {
      expected_files <- vapply(res, function(item) item$result_file, character(1))
      if (any(!file.exists(expected_files))) {
        stop("Successful workers did not write the expected result files: ",
             paste(expected_files[!file.exists(expected_files)], collapse = ", "))
      }
    }

    merge_info <- NULL
    if (isTRUE(auto_merge)) {
      msg("Starting automatic merge of per-unit results...")
      merge_info <- merge_dss_results(
        results_dir = results_dir, result_prefix = result_prefix, merged_dir = outdir,
        p_col = "pvals", old_fdr_col = "fdrs", new_fdr_col = "fdr_global",
        overwrite_fdr = TRUE, save_merged_rds = TRUE, save_merged_tsv = TRUE
      )
      msg("Merged results written to: ", merge_info$merged_file)
    }
    
    invisible(list(
      mode = "live", summary = res_df, summary_file = summary_file,
      manifest_file = manifest_file, merge_info = merge_info, outdir = outdir
    ))
  }
  
  # ----------------------------
  # slurm backend writer
  # ----------------------------
  write_slurm_backend <- function() {
    msg("Writing Slurm backend assets to: ", outdir)
    
    # Save metadata once.
    safe_save_rds(meta2, meta_file)
    
    # Pre-write a SMALL subset BSseq object per work unit, and a slurm
    # manifest that maps task_id -> label, chr, unit_file. Each array task
    # then loads ONLY its block (huge memory reduction vs full object).
    msg("Writing ", length(work_units), " per-unit subset objects to: ", units_dir)
    unit_files <- character(length(work_units))
    for (i in seq_along(work_units)) {
      u <- work_units[[i]]
      uf <- file.path(units_dir, paste0(result_prefix, "_", u$label, "_bs.rds"))
      saveRDS(bs[u$idx, ], uf, compress = slurm_subset_compress)
      unit_files[i] <- uf
    }
    
    slurm_manifest <- data.frame(
      task_id   = seq_along(work_units),
      label     = vapply(work_units, function(u) u$label, character(1)),
      chr       = vapply(work_units, function(u) u$chr,   character(1)),
      n_loci    = vapply(work_units, function(u) length(u$idx), integer(1)),
      unit_file = unit_files,
      stringsAsFactors = FALSE
    )
    slurm_manifest_file <- file.path(meta_dir, paste0(result_prefix, "_slurm_manifest.tsv"))
    write.table(slurm_manifest, slurm_manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    worker_script <- file.path(scripts_dir, paste0(result_prefix, "_worker.R"))
    sbatch_script <- file.path(scripts_dir, paste0(result_prefix, "_array.sbatch"))
    merge_script  <- file.path(scripts_dir, paste0(result_prefix, "_merge.R"))
    merge_sbatch  <- file.path(scripts_dir, paste0(result_prefix, "_merge.sbatch"))
    
    # serialize contrast structures if present
    contrast_file <- file.path(meta_dir, paste0(result_prefix, "_contrast_structure.rds"))
    if (!is.null(contrast_structure)) safe_save_rds(contrast_structure, contrast_file) else contrast_file <- NA_character_
    manual_contrast_file <- file.path(meta_dir, paste0(result_prefix, "_manual_contrast.rds"))
    if (!is.null(Contrast)) safe_save_rds(Contrast, manual_contrast_file) else manual_contrast_file <- NA_character_
    
    # ---- worker.R : loads ONE pre-subset block, fits, tests, saves ----
    worker_code <- c(
      '## auto-generated DSS block worker',
      'Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",',
      '           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")',
      'if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {',
      '  try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)',
      '  try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)',
      '}',
      'args <- commandArgs(trailingOnly = TRUE)',
      'arg_list <- list()',
      'for (a in args) {',
      '  sp <- strsplit(a, "=", fixed = TRUE)[[1]]',
      '  key <- sp[1]; val <- paste(sp[-1], collapse = "=")',
      '  arg_list[[key]] <- val',
      '}',
      '`%||%` <- function(x, y) if (is.null(x) || identical(x, "NA")) y else x',
      '',
      'load_pkg <- function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)',
      'for (pkg in c("DSS","bsseq","GenomicRanges","SummarizedExperiment")) load_pkg(pkg)',
      '',
      'normalize_char_levels <- function(x) trimws(as.character(x))',
      '',
      'validate_contrast_structure <- function(cs, levels_outcome) {',
      '  if (!is.list(cs) || is.null(names(cs)) || any(names(cs) == "")) stop("Bad contrast_structure")',
      '  for (nm in names(cs)) {',
      '    one <- cs[[nm]]',
      '    if (!is.list(one) || !all(c("case","control") %in% names(one))) stop("Bad contrast entry: ", nm)',
      '    case <- normalize_char_levels(one$case); control <- normalize_char_levels(one$control)',
      '    if (any(!case %in% levels_outcome)) stop("Unknown case levels in ", nm)',
      '    if (any(!control %in% levels_outcome)) stop("Unknown control levels in ", nm)',
      '    if (length(intersect(case, control)) > 0) stop("Overlap in contrast ", nm)',
      '  }',
      '  invisible(TRUE)',
      '}',
      '',
      'coerce_outcome_type_auto <- function(x) {',
      '  x_non_na <- x[!is.na(x)]',
      '  if (is.numeric(x_non_na) || is.integer(x_non_na)) return("continuous")',
      '  if (is.logical(x_non_na)) return("binary")',
      '  if (is.factor(x_non_na) || is.character(x_non_na)) {',
      '    u <- unique(as.character(x_non_na))',
      '    if (length(u) == 2) return("binary")',
      '    return("categorical")',
      '  }',
      '  "categorical"',
      '}',
      '',
      'get_memory_snapshot_mb <- function() {',
      '  rss_mb <- NA_real_; vmhwm_mb <- NA_real_',
      '  status_file <- "/proc/self/status"',
      '  if (file.exists(status_file)) {',
      '    txt <- tryCatch(readLines(status_file, warn = FALSE), error = function(e) character())',
      '    if (length(txt)) {',
      '      rss_line <- grep("^VmRSS:", txt, value = TRUE)',
      '      hwm_line <- grep("^VmHWM:", txt, value = TRUE)',
      '      parse_kb <- function(z) { if (!length(z)) return(NA_real_); as.numeric(gsub("[^0-9]", "", z)) / 1024 }',
      '      rss_mb <- parse_kb(rss_line); vmhwm_mb <- parse_kb(hwm_line)',
      '    }',
      '  }',
      '  gc_mat <- gc(); gc_used_mb <- sum(gc_mat[, "used"] * c(8, 8)) / 1024^2',
      '  data.frame(rss_mb = rss_mb, vmhwm_mb = vmhwm_mb, gc_used_mb = gc_used_mb)',
      '}',
      '',
      'write_log_line <- function(log_file, ...) {',
      '  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))',
      '  cat(line, "\\n", file = log_file, append = TRUE)',
      '}',
      '',
      'build_test_plan <- function(meta_df, formula, outcome_col, outcome_type_final, ref,',
      '                            contrast_structure, multi_group_default, coef, term, Contrast) {',
      '  n_manual <- sum(!is.null(coef), !is.null(term), !is.null(Contrast))',
      '  manual_mode <- n_manual > 0',
      '  if (manual_mode) {',
      '    return(list(mode = "manual", test_plan = list(list(',
      '      name = "manual",',
      '      type = if (!is.null(coef)) "coef" else if (!is.null(term)) "term" else "contrast",',
      '      coef = coef, term = term, Contrast = Contrast',
      '    )), meta = meta_df, formula = formula, outcome_type_final = outcome_type_final))',
      '  }',
      '  term_labels <- attr(stats::terms(formula), "term.labels")',
      '  cov_terms <- setdiff(term_labels, outcome_col)',
      '  if (outcome_type_final == "continuous") {',
      '    return(list(mode = "auto", test_plan = list(list(name = "main", type = "coef", coef = outcome_col)),',
      '                meta = meta_df, formula = stats::reformulate(c(outcome_col, cov_terms)),',
      '                outcome_type_final = outcome_type_final))',
      '  }',
      '  if (outcome_type_final == "binary") {',
      '    x <- meta_df[[outcome_col]]',
      '    if (is.logical(x)) x <- ifelse(x, "TRUE", "FALSE")',
      '    x <- factor(normalize_char_levels(x)); lv <- levels(x)',
      '    if (length(lv) != 2) stop("Binary outcome requires 2 levels.")',
      '    if (is.null(ref) || !ref %in% lv) stop("Invalid ref.")',
      '    x <- stats::relevel(x, ref = ref); meta_df[[outcome_col]] <- x',
      '    return(list(mode = "auto", test_plan = list(list(name = "main", type = "coef", coef = outcome_col)),',
      '                meta = meta_df, formula = stats::reformulate(c(outcome_col, cov_terms)),',
      '                outcome_type_final = outcome_type_final))',
      '  }',
      '  if (outcome_type_final == "categorical") {',
      '    x <- factor(normalize_char_levels(meta_df[[outcome_col]])); lv <- levels(x)',
      '    if (length(lv) < 3) stop("Categorical outcome requires 3+ levels.")',
      '    meta_df[[outcome_col]] <- x',
      '    formula_fit <- stats::reformulate(c(outcome_col, cov_terms)); plan <- list()',
      '    if (multi_group_default == "term") {',
      '      plan[[length(plan) + 1L]] <- list(name = paste0(outcome_col, "_omnibus"), type = "term", term = outcome_col)',
      '    }',
      '    if (multi_group_default == "one_vs_rest") {',
      '      if (is.null(contrast_structure)) {',
      '        contrast_structure <- setNames(lapply(lv, function(L) list(case = L, control = setdiff(lv, L))),',
      '                                      paste0(lv, "_vs_rest"))',
      '      }',
      '    }',
      '    if (!is.null(contrast_structure)) {',
      '      validate_contrast_structure(contrast_structure, lv)',
      '      for (nm in names(contrast_structure)) {',
      '        plan[[length(plan) + 1L]] <- list(name = nm, type = "level_contrast",',
      '          case = normalize_char_levels(contrast_structure[[nm]]$case),',
      '          control = normalize_char_levels(contrast_structure[[nm]]$control), outcome_col = outcome_col)',
      '      }',
      '    }',
      '    if (!length(plan)) stop("No categorical test plan created.")',
      '    return(list(mode = "auto", test_plan = plan, meta = meta_df, formula = formula_fit,',
      '                outcome_type_final = outcome_type_final))',
      '  }',
      '  stop("Unhandled outcome type")',
      '}',
      '',
      'build_factor_level_contrast <- function(fit, outcome_col, case, control, meta_df) {',
      '  xcols <- colnames(fit$X); levs <- levels(meta_df[[outcome_col]]); ref_level <- levs[1]',
      '  coef_for_level <- function(level_name) {',
      '    if (identical(level_name, ref_level)) return(NULL)',
      '    exact1 <- paste0(outcome_col, level_name); exact2 <- paste0(outcome_col, make.names(level_name))',
      '    hits <- which(xcols %in% c(exact1, exact2)); if (length(hits) == 1) return(hits)',
      '    hits <- grep(paste0("^", gsub("([.|()\\\\^{}+$*?]|\\\\[|\\\\])", "\\\\\\\\\\1", outcome_col)), xcols)',
      '    hit_names <- xcols[hits]',
      '    hit_idx <- hits[grepl(make.names(level_name), make.names(hit_names), fixed = TRUE) |',
      '                    grepl(level_name, hit_names, fixed = TRUE)]',
      '    if (length(hit_idx) == 1) return(hit_idx)',
      '    stop("Could not map level ", level_name, " to fit$X columns.")',
      '  }',
      '  v <- rep(0, ncol(fit$X)); w_case <- 1 / length(case); w_ctrl <- 1 / length(control)',
      '  for (lv in case) { idx <- coef_for_level(lv); if (!is.null(idx)) v[idx] <- v[idx] + w_case }',
      '  for (lv in control) { idx <- coef_for_level(lv); if (!is.null(idx)) v[idx] <- v[idx] - w_ctrl }',
      '  matrix(v, ncol = 1)',
      '}',
      '',
      '## ---- resolve this task from the manifest ----',
      'manifest_file <- arg_list[["manifest_file"]]',
      'task_id <- as.integer(arg_list[["task_id"]])',
      'mani <- utils::read.table(manifest_file, header = TRUE, sep = "\\t", stringsAsFactors = FALSE)',
      'row <- mani[mani$task_id == task_id, , drop = FALSE]',
      'if (nrow(row) != 1) stop("task_id ", task_id, " not found uniquely in manifest")',
      'label <- row$label; chr <- row$chr; unit_file <- row$unit_file',
      '',
      'meta_file <- arg_list[["meta_file"]]',
      'formula_text <- arg_list[["formula_text"]]',
      'outcome_col <- arg_list[["outcome_col"]]',
      'outcome_type <- arg_list[["outcome_type"]]',
      'ref <- if (identical(arg_list[["ref"]], "NA")) NULL else arg_list[["ref"]]',
      'coef <- if (identical(arg_list[["coef"]], "NA")) NULL else arg_list[["coef"]]',
      'term <- if (identical(arg_list[["term"]], "NA")) NULL else arg_list[["term"]]',
      'contrast_file <- if (identical(arg_list[["contrast_file"]], "NA")) NULL else arg_list[["contrast_file"]]',
      'manual_contrast_file <- if (identical(arg_list[["manual_contrast_file"]], "NA")) NULL else arg_list[["manual_contrast_file"]]',
      'multi_group_default <- arg_list[["multi_group_default"]]',
      'smoothing <- as.logical(arg_list[["smoothing"]])',
      'smoothing.span <- as.numeric(arg_list[["smoothing.span"]])',
      'results_dir <- arg_list[["results_dir"]]',
      'logs_dir <- arg_list[["logs_dir"]]',
      'result_prefix <- arg_list[["result_prefix"]]',
      'save_fit <- as.logical(arg_list[["save_fit"]])',
      '',
      'log_file <- file.path(logs_dir, paste0(result_prefix, "_", label, ".log"))',
      'write_log_line(log_file, "Starting worker for ", label, " (chr=", chr, ")")',
      '',
      'bs_chr <- readRDS(unit_file)   ## already a small per-block subset',
      'meta <- readRDS(meta_file)',
      'formula <- stats::as.formula(formula_text)',
      'contrast_structure <- if (!is.null(contrast_file)) readRDS(contrast_file) else NULL',
      'Contrast <- if (!is.null(manual_contrast_file)) readRDS(manual_contrast_file) else NULL',
      'outcome_type_final <- if (identical(outcome_type, "auto")) coerce_outcome_type_auto(meta[[outcome_col]]) else outcome_type',
      '',
      'n_loci <- nrow(bs_chr)',
      'mem0 <- get_memory_snapshot_mb()',
      'write_log_line(log_file, "Initial memory MB: rss=", round(mem0$rss_mb,1), ", vmhwm=", round(mem0$vmhwm_mb,1))',
      'write_log_line(log_file, "n loci: ", n_loci)',
      '',
      'test_plan_obj <- build_test_plan(meta_df = meta, formula = formula, outcome_col = outcome_col,',
      '  outcome_type_final = outcome_type_final, ref = ref, contrast_structure = contrast_structure,',
      '  multi_group_default = multi_group_default, coef = coef, term = term, Contrast = Contrast)',
      'meta_use <- test_plan_obj$meta; formula_use <- test_plan_obj$formula',
      'X_preview <- model.matrix(formula_use, data = meta_use)',
      'if (qr(X_preview)$rank < ncol(X_preview)) stop("Rank-deficient design on ", label)',
      '',
      't0_fit <- proc.time()[["elapsed"]]',
      'fit <- DSS::DMLfit.multiFactor(BSobj = bs_chr, design = meta_use, formula = formula_use,',
      '                               smoothing = smoothing, smoothing.span = smoothing.span)',
      't_fit <- proc.time()[["elapsed"]] - t0_fit',
      'write_log_line(log_file, "Fit done in sec: ", round(t_fit,2))',
      'write_log_line(log_file, "fit$X cols: ", paste(colnames(fit$X), collapse = ", "))',
      '',
      'tests_out <- list()',
      'for (tp in test_plan_obj$test_plan) {',
      '  write_log_line(log_file, "Running test: ", tp$name, " [", tp$type, "]")',
      '  t0 <- proc.time()[["elapsed"]]',
      '  test_res <- switch(tp$type,',
      '    coef = DSS::DMLtest.multiFactor(fit, coef = tp$coef),',
      '    term = DSS::DMLtest.multiFactor(fit, term = tp$term),',
      '    contrast = DSS::DMLtest.multiFactor(fit, Contrast = tp$Contrast),',
      '    level_contrast = {',
      '      Cmat <- build_factor_level_contrast(fit, tp$outcome_col, tp$case, tp$control, meta_use)',
      '      DSS::DMLtest.multiFactor(fit, Contrast = Cmat)',
      '    },',
      '    stop("Unknown test type: ", tp$type))',
      '  dt <- proc.time()[["elapsed"]] - t0',
      '  write_log_line(log_file, "Test done: ", tp$name, " sec=", round(dt,2), " loci=", nrow(test_res))',
      '  tests_out[[tp$name]] <- test_res   ## bare table -> same structure as live backend',
      '}',
      '',
      'mem1 <- get_memory_snapshot_mb()',
      'res_file <- file.path(results_dir, paste0(result_prefix, "_", label, ".rds"))',
      'out <- list(chr = chr, label = label, status = "ok", n_loci = n_loci,',
      '            tests = tests_out, fit = if (save_fit) fit else NULL,',
      '            timing_sec = list(fit = t_fit), memory_mb = list(start = mem0, end = mem1))',
      'saveRDS(out, res_file, compress = "gzip")',
      'write_log_line(log_file, "Saved: ", res_file)',
      'write_log_line(log_file, "Final memory MB: rss=", round(mem1$rss_mb,1), ", vmhwm=", round(mem1$vmhwm_mb,1))',
      'write_log_line(log_file, "Worker complete for ", label)'
    )
    writeLines(worker_code, con = worker_script)
    
    # ---- merge.R : reproduces the live merged output (global BH FDR) ----
    merge_code <- c(
      '## auto-generated merge: genome-wide BH FDR across all blocks',
      'args <- commandArgs(trailingOnly = TRUE)',
      'arg_list <- list()',
      'for (a in args) { sp <- strsplit(a, "=", fixed = TRUE)[[1]]; arg_list[[sp[1]]] <- paste(sp[-1], collapse = "=") }',
      'results_dir <- arg_list[["results_dir"]]',
      'merged_dir  <- arg_list[["merged_dir"]]',
      'result_prefix <- arg_list[["result_prefix"]]',
      '',
      'files <- list.files(results_dir, pattern = paste0("^", result_prefix, "_.*\\\\.rds$"), full.names = TRUE)',
      'if (!length(files)) stop("No per-unit result files found.")',
      'res_ok <- Filter(function(x) identical(x$status, "ok"), lapply(files, readRDS))',
      'if (!length(res_ok)) stop("No successful result files found.")',
      'shared <- Reduce(intersect, lapply(res_ok, function(x) names(x$tests)))',
      'if (!length(shared)) stop("No shared test names across files.")',
      '',
      'merged_tables <- list()',
      'for (nm in shared) {',
      '  tabs <- lapply(res_ok, function(x) { tab <- x$tests[[nm]]; if (is.null(tab)) return(NULL); tab$chr_job <- x$chr; tab })',
      '  tab <- do.call(rbind, Filter(Negate(is.null), tabs)); rownames(tab) <- NULL',
      '  if ("fdrs" %in% colnames(tab)) tab$fdrs_per_chr <- tab$fdrs',
      '  tab$fdr_global <- stats::p.adjust(tab$pvals, method = "BH")',
      '  tab$fdrs <- tab$fdr_global',
      '  merged_tables[[nm]] <- tab',
      '  con <- gzfile(file.path(merged_dir, paste0(result_prefix, "_", nm, "_merged.tsv.gz")), "wt")',
      '  utils::write.table(tab, con, sep = "\\t", row.names = FALSE, quote = FALSE); close(con)',
      '}',
      'saveRDS(merged_tables, file.path(merged_dir, paste0(result_prefix, "_merged_results.rds")), compress = "gzip")',
      'summary_df <- do.call(rbind, lapply(res_ok, function(x) data.frame(chr = x$chr, label = x$label, n_loci = x$n_loci, stringsAsFactors = FALSE)))',
      'utils::write.table(summary_df, file.path(merged_dir, paste0(result_prefix, "_merged_summary.tsv")), sep = "\\t", row.names = FALSE, quote = FALSE)',
      'cat("Merge complete:", file.path(merged_dir, paste0(result_prefix, "_merged_results.rds")), "\\n")'
    )
    writeLines(merge_code, con = merge_script)
    
    formula_text <- paste(deparse(formula), collapse = "")
    array_n <- nrow(slurm_manifest)
    
    sbatch_lines <- c(
      "#!/bin/bash",
      paste0("#SBATCH --job-name=", result_prefix),
      if (!is.null(slurm_partition)) paste0("#SBATCH --partition=", slurm_partition),
      paste0("#SBATCH --cpus-per-task=", slurm_cpus_per_task),
      paste0("#SBATCH --mem=", slurm_mem),
      paste0("#SBATCH --time=", slurm_time),
      paste0("#SBATCH --array=1-", array_n),
      paste0("#SBATCH --output=", file.path(logs_dir, paste0(result_prefix, "_slurm_%A_%a.out"))),
      paste0("#SBATCH --error=", file.path(logs_dir, paste0(result_prefix, "_slurm_%A_%a.err")))
    )
    if (!is.null(slurm_account) && nzchar(slurm_account)) {
      sbatch_lines <- c(sbatch_lines, paste0("#SBATCH --account=", slurm_account))
    }
    sbatch_lines <- c(
      sbatch_lines, "",
      "set -euo pipefail", "",
      "# pin BLAS/OpenMP so each task uses exactly one core",
      "export OMP_NUM_THREADS=1",
      "export OPENBLAS_NUM_THREADS=1",
      "export MKL_NUM_THREADS=1",
      "export VECLIB_MAXIMUM_THREADS=1", "",
      slurm_module_load, "",
      paste0("MANIFEST='", slurm_manifest_file, "'"),
      paste0("/usr/bin/time -v ", rscript_path, " ", worker_script, " \\"),
      paste0("  manifest_file='", slurm_manifest_file, "' \\"),
      paste0("  task_id=$SLURM_ARRAY_TASK_ID \\"),
      paste0("  meta_file='", meta_file, "' \\"),
      paste0("  formula_text='", formula_text, "' \\"),
      paste0("  outcome_col='", outcome_col, "' \\"),
      paste0("  outcome_type='", outcome_type, "' \\"),
      paste0("  ref='", ifelse(is.null(ref), "NA", ref), "' \\"),
      paste0("  coef='", ifelse(is.null(coef), "NA", as.character(coef)), "' \\"),
      paste0("  term='", ifelse(is.null(term), "NA", term), "' \\"),
      paste0("  contrast_file='", ifelse(is.na(contrast_file), "NA", contrast_file), "' \\"),
      paste0("  manual_contrast_file='", ifelse(is.na(manual_contrast_file), "NA", manual_contrast_file), "' \\"),
      paste0("  multi_group_default='", multi_group_default, "' \\"),
      paste0("  smoothing='", ifelse(smoothing, "TRUE", "FALSE"), "' \\"),
      paste0("  smoothing.span='", smoothing.span, "' \\"),
      paste0("  results_dir='", results_dir, "' \\"),
      paste0("  logs_dir='", logs_dir, "' \\"),
      paste0("  result_prefix='", result_prefix, "' \\"),
      paste0("  save_fit='", ifelse(save_fit, "TRUE", "FALSE"), "'")
    )
    writeLines(sbatch_lines, con = sbatch_script)
    
    # ---- merge.sbatch : single task, runs after the array succeeds ----
    merge_sbatch_lines <- c(
      "#!/bin/bash",
      paste0("#SBATCH --job-name=", result_prefix, "_merge"),
      if (!is.null(slurm_partition)) paste0("#SBATCH --partition=", slurm_partition),
      "#SBATCH --cpus-per-task=1",
      paste0("#SBATCH --mem=", slurm_mem),
      paste0("#SBATCH --time=", slurm_time),
      paste0("#SBATCH --output=", file.path(logs_dir, paste0(result_prefix, "_merge_%j.out"))),
      paste0("#SBATCH --error=",  file.path(logs_dir, paste0(result_prefix, "_merge_%j.err")))
    )
    if (!is.null(slurm_account) && nzchar(slurm_account)) {
      merge_sbatch_lines <- c(merge_sbatch_lines, paste0("#SBATCH --account=", slurm_account))
    }
    merge_sbatch_lines <- c(
      merge_sbatch_lines, "",
      "set -euo pipefail", "",
      slurm_module_load, "",
      paste0(rscript_path, " ", merge_script, " \\"),
      paste0("  results_dir='", results_dir, "' \\"),
      paste0("  merged_dir='", outdir, "' \\"),
      paste0("  result_prefix='", result_prefix, "'")
    )
    writeLines(merge_sbatch_lines, con = merge_sbatch)
    
    Sys.chmod(worker_script, mode = "0755")
    Sys.chmod(sbatch_script, mode = "0755")
    Sys.chmod(merge_script, mode = "0755")
    Sys.chmod(merge_sbatch, mode = "0755")
    
    instr_file <- file.path(outdir, paste0(result_prefix, "_slurm_instructions.txt"))
    writeLines(c(
      "# Submit the array, then the merge (dependency: merge runs only if array succeeds):",
      paste0("AID=$(sbatch --parsable ", sbatch_script, ")"),
      paste0("sbatch --dependency=afterok:${AID} ", merge_sbatch),
      "",
      paste0("Array size: ", array_n, " tasks"),
      paste0("Slurm manifest: ", slurm_manifest_file),
      paste0("Per-block subset objects: ", units_dir),
      paste0("Worker script: ", worker_script),
      paste0("Merge script: ", merge_script),
      paste0("Meta RDS: ", meta_file),
      paste0("Results dir: ", results_dir),
      paste0("Merged output dir: ", outdir),
      paste0("Logs dir: ", logs_dir)
    ), con = instr_file)
    
    msg("Wrote worker: ", worker_script)
    msg("Wrote array sbatch: ", sbatch_script)
    msg("Wrote merge script + sbatch: ", merge_script)
    msg("Array size: ", array_n, " tasks. Submit instructions: ", instr_file)
    
    invisible(list(
      mode = "slurm",
      worker_script = worker_script,
      sbatch_script = sbatch_script,
      merge_script = merge_script,
      merge_sbatch = merge_sbatch,
      slurm_manifest_file = slurm_manifest_file,
      units_dir = units_dir,
      meta_file = meta_file,
      results_dir = results_dir,
      logs_dir = logs_dir,
      instructions_file = instr_file,
      outdir = outdir
    ))
  }
  
  # ----------------------------
  # dispatch
  # ----------------------------
  if (backend == "live") {
    return(run_live_backend())
  } else {
    return(write_slurm_backend())
  }
}
