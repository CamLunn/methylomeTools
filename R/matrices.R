
# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
.to_numeric_matrix <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (!is.matrix(x)) x <- as.matrix(x)
  suppressWarnings(mode(x) <- "numeric")
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(name, " must have rownames (loci) and colnames (samples).")
  }
  x
}

.align_cov_and_mats <- function(cov_mat, mats) {
  mats <- mats[!vapply(mats, is.null, logical(1))]
  if (length(mats) == 0) {
    stop("At least one matrix besides cov_mat must be supplied.")
  }
  
  loci_common <- rownames(cov_mat)
  samp_common <- colnames(cov_mat)
  
  for (nm in names(mats)) {
    loci_common <- intersect(loci_common, rownames(mats[[nm]]))
    samp_common <- intersect(samp_common, colnames(mats[[nm]]))
  }
  
  if (length(loci_common) == 0) stop("No overlapping loci across supplied matrices.")
  if (length(samp_common) == 0) stop("No overlapping samples across supplied matrices.")
  
  list(
    cov_mat = cov_mat[loci_common, samp_common, drop = FALSE],
    mats = lapply(mats, function(m) m[loci_common, samp_common, drop = FALSE])
  )
}

.compute_beta_from_meth <- function(meth_mat, cov_mat) {
  beta <- meth_mat / cov_mat
  beta[cov_mat <= 0 | is.na(cov_mat)] <- NA_real_
  beta
}

.make_col_quantile_matrix_filtered <- function(cov_mat,
                                               probs,
                                               lo.count = NULL,
                                               hi.count = NULL) {
  qv <- vapply(
    seq_len(ncol(cov_mat)),
    function(j) {
      x <- cov_mat[, j]
      x <- x[!is.na(x) & is.finite(x) & x > 0]
      
      if (!is.null(lo.count)) {
        x <- x[x >= lo.count]
      }
      
      if (!is.null(hi.count)) {
        x <- x[x <= hi.count]
      }
      
      if (!length(x)) {
        return(NA_real_)
      }
      
      stats::quantile(
        x,
        probs = probs,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    },
    numeric(1)
  )
  
  matrix(
    rep(qv, each = nrow(cov_mat)),
    nrow = nrow(cov_mat),
    ncol = ncol(cov_mat),
    dimnames = dimnames(cov_mat)
  )
}

.apply_row_subset_rule <- function(pass_mat,
                                   keep_rule = c("any", "all", "max_missing"),
                                   max_missing_prop = NULL) {
  keep_rule <- match.arg(keep_rule)
  
  switch(
    keep_rule,
    any = rowSums(pass_mat, na.rm = TRUE) >= 1,
    all = rowSums(pass_mat, na.rm = TRUE) == ncol(pass_mat),
    max_missing = {
      if (is.null(max_missing_prop)) {
        rep(TRUE, nrow(pass_mat))
      } else {
        rowMeans(!pass_mat, na.rm = TRUE) <= max_missing_prop
      }
    }
  )
}

.format_bin_name <- function(lower, upper, prefix = "cov") {
  if (is.null(upper) || !is.finite(upper)) {
    paste0(prefix, lower, "_Inf")
  } else {
    paste0(prefix, lower, "_", upper)
  }
}

# ------------------------------------------------------------
# Build one combined pass matrix from all requested filters
# ------------------------------------------------------------
.build_combined_pass_matrix <- function(cov_mat,
                                        lo.count = NULL,
                                        hi.count = NULL,
                                        lo.perc = NULL,
                                        hi.perc = NULL,
                                        percentile_after_count_filters = TRUE) {
  pass <- matrix(TRUE, nrow(cov_mat), ncol(cov_mat), dimnames = dimnames(cov_mat))
  
  if (!is.null(lo.count)) {
    pass <- pass & (cov_mat >= lo.count)
  }
  
  if (!is.null(hi.count)) {
    pass <- pass & (cov_mat <= hi.count)
  }
  
  q_lo_count <- if (percentile_after_count_filters) lo.count else NULL
  q_hi_count <- if (percentile_after_count_filters) hi.count else NULL
  
  if (!is.null(lo.perc)) {
    lo_mat <- .make_col_quantile_matrix_filtered(
      cov_mat = cov_mat,
      probs = lo.perc / 100,
      lo.count = q_lo_count,
      hi.count = q_hi_count
    )
    pass <- pass & (cov_mat >= lo_mat)
  }
  
  if (!is.null(hi.perc)) {
    hi_mat <- .make_col_quantile_matrix_filtered(
      cov_mat = cov_mat,
      probs = hi.perc / 100,
      lo.count = q_lo_count,
      hi.count = q_hi_count
    )
    pass <- pass & (cov_mat <= hi_mat)
  }
  
  pass[is.na(cov_mat) | !is.finite(cov_mat)] <- FALSE
  pass
}

# ------------------------------------------------------------
# Generic one-shot matrix masking/subsetting function
# ------------------------------------------------------------
mask_matrices_by_coverage <- function(
    cov_mat,
    beta_mat = NULL,
    meth_mat = NULL,
    extra_mats = NULL,
    lo.count = NULL,
    hi.count = NULL,
    lo.perc = NULL,
    hi.perc = NULL,
    mode = c("mask", "subset"),
    keep_rule = c("any", "all", "max_missing"),
    max_missing_prop = 1.0,
    return_what = c("beta", "meth", "cov", "all"),
    verbose = TRUE
) {
  mode <- match.arg(mode)
  keep_rule <- match.arg(keep_rule)
  return_what <- match.arg(return_what)
  
  cov_mat  <- .to_numeric_matrix(cov_mat, "cov_mat")
  beta_mat <- .to_numeric_matrix(beta_mat, "beta_mat")
  meth_mat <- .to_numeric_matrix(meth_mat, "meth_mat")
  
  if (is.null(beta_mat) && !is.null(meth_mat)) {
    beta_mat <- .compute_beta_from_meth(meth_mat, cov_mat)
  }
  
  mats <- c(
    list(beta = beta_mat, meth = meth_mat, cov = cov_mat),
    extra_mats
  )
  
  aligned <- .align_cov_and_mats(cov_mat = cov_mat, mats = mats)
  cov  <- aligned$cov_mat
  mats <- aligned$mats
  
  pass_mat <- .build_combined_pass_matrix(
    cov_mat  = cov,
    lo.count = lo.count,
    hi.count = hi.count,
    lo.perc  = lo.perc,
    hi.perc  = hi.perc
  )
  
  if (verbose) {
    n_entries_total <- length(pass_mat)
    n_entries_keep  <- sum(pass_mat, na.rm = TRUE)
    pct_entries <- round(100 * n_entries_keep / n_entries_total, 2)
    
    n_cpg_any <- sum(rowSums(pass_mat, na.rm = TRUE) > 0)
    pct_cpg_any <- round(100 * n_cpg_any / nrow(pass_mat), 2)
    
    message(
      "Coverage filter retained ", n_entries_keep, " / ", n_entries_total,
      " matrix entries (", pct_entries, "%); ",
      n_cpg_any, " / ", nrow(pass_mat), " CpGs have at least one passing sample (",
      pct_cpg_any, "%)."
    )
  }
  
  if (mode == "mask") {
    out <- lapply(mats, function(m) {
      m2 <- m
      m2[!pass_mat] <- NA_real_
      m2
    })
    
    if (verbose) {
      n_cpg_nonallna <- sum(rowSums(!is.na(out$cov)) > 0)
      message(
        "Mask mode output retains ", n_cpg_nonallna, " / ", nrow(out$cov),
        " CpGs with at least one non-NA value."
      )
    }
  } else {
    keep_rows <- .apply_row_subset_rule(
      pass_mat = pass_mat,
      keep_rule = keep_rule,
      max_missing_prop = max_missing_prop
    )
    
    out <- lapply(mats, function(m) m[keep_rows, , drop = FALSE])
    
    if (verbose) {
      pct_rows <- round(100 * sum(keep_rows) / length(keep_rows), 2)
      message(
        "Subset mode retained ", sum(keep_rows), " / ", length(keep_rows),
        " CpGs (", pct_rows, "%) using keep_rule='", keep_rule, "'."
      )
    }
  }
  
  if (return_what == "all") return(out)
  out[[return_what]]
}


# ------------------------------------------------------------
# Threshold series wrapper:
# - default: cumulative lower bounds [thr, Inf)
# - sections: interval bins [thr_i, thr_{i+1})
# All thresholds can also be combined with lo/hi count/perc filters
# ------------------------------------------------------------
make_coverage_threshold_matrices <- function(
    beta_base = NULL,
    meth_mat = NULL,
    cov_mat,
    thresholds = c(5, 10, 20, 30, 50, 100),
    cut_threshold = c("default", "sections"),
    mode = c("subset", "mask"),
    keep_rule = c("any", "all", "max_missing"),
    max_missing_prop = NULL,
    lo.count = NULL,
    hi.count = NULL,
    lo.perc = NULL,
    hi.perc = NULL,
    return_what = c("beta", "meth", "cov", "all"),
    extra_mats = NULL,
    verbose = TRUE
) {
  cut_threshold <- match.arg(cut_threshold)
  mode <- match.arg(mode)
  keep_rule <- match.arg(keep_rule)
  return_what <- match.arg(return_what)
  
  cov_mat   <- .to_numeric_matrix(cov_mat, "cov_mat")
  beta_base <- .to_numeric_matrix(beta_base, "beta_base")
  meth_mat  <- .to_numeric_matrix(meth_mat, "meth_mat")
  
  if (is.null(beta_base) && is.null(meth_mat) && return_what %in% c("beta", "all")) {
    stop("Supply either beta_base or meth_mat.")
  }
  
  thresholds <- sort(unique(as.numeric(thresholds)))
  if (length(thresholds) == 0 || anyNA(thresholds) || any(thresholds < 0)) {
    stop("thresholds must be non-negative numeric values.")
  }
  
  if (!is.null(lo.perc) && (lo.perc < 0 || lo.perc > 100)) stop("lo.perc must be in [0,100].")
  if (!is.null(hi.perc) && (hi.perc < 0 || hi.perc > 100)) stop("hi.perc must be in [0,100].")
  if (!is.null(lo.perc) && !is.null(hi.perc) && lo.perc > hi.perc) {
    stop("lo.perc cannot be greater than hi.perc.")
  }
  if (!is.null(lo.count) && !is.null(hi.count) && lo.count > hi.count) {
    stop("lo.count cannot be greater than hi.count.")
  }
  if (!is.null(max_missing_prop) &&
      (!is.numeric(max_missing_prop) || length(max_missing_prop) != 1L ||
       is.na(max_missing_prop) || max_missing_prop < 0 || max_missing_prop > 1)) {
    stop("max_missing_prop must be a single numeric value in [0, 1].")
  }
  
  bins <- if (cut_threshold == "default") {
    lapply(thresholds, function(thr) list(lower = thr, upper = Inf))
  } else {
    uppers <- c(thresholds[-1], Inf)
    Map(function(lo, hi) list(lower = lo, upper = hi), thresholds, uppers)
  }
  
  out <- vector("list", length(bins))
  names(out) <- vapply(
    bins,
    function(b) .format_bin_name(b$lower, b$upper, prefix = "mincov"),
    character(1)
  )
  
  for (i in seq_along(bins)) {
    b <- bins[[i]]
    lo_count_i <- b$lower
    hi_count_i <- if (is.finite(b$upper)) b$upper else NULL
    
    if (verbose) {
      rng <- if (is.null(hi_count_i)) {
        paste0("[", lo_count_i, ", Inf)")
      } else {
        paste0("[", lo_count_i, ", ", hi_count_i, ")")
      }
      message("Building ", names(out)[i], " using coverage bin ", rng)
      msg_parts <- c(
        paste0("  base bin: ", rng),
        if (!is.null(lo.count)) paste0("extra lo.count=", lo.count) else NULL,
        if (!is.null(hi.count)) paste0("extra hi.count=", hi.count) else NULL,
        if (!is.null(lo.perc))  paste0("extra lo.perc=", lo.perc)   else NULL,
        if (!is.null(hi.perc))  paste0("extra hi.perc=", hi.perc)   else NULL,
        if (!is.null(max_missing_prop)) paste0("max_missing_prop=", max_missing_prop) else NULL
      )
      message(paste(msg_parts, collapse = " | "))
    }
    
    eff_lo_count <- lo_count_i
    if (!is.null(lo.count)) eff_lo_count <- max(eff_lo_count, lo.count)
    
    eff_hi_count <- hi_count_i
    if (!is.null(hi.count)) {
      eff_hi_count <- if (is.null(eff_hi_count)) hi.count else min(eff_hi_count, hi.count)
    }
    
    if (!is.null(eff_hi_count) && cut_threshold == "sections") {
      eff_hi_count_internal <- eff_hi_count - .Machine$double.eps^0.5
    } else {
      eff_hi_count_internal <- eff_hi_count
    }
    
    out[[i]] <- mask_matrices_by_coverage(
      cov_mat = cov_mat,
      beta_mat = beta_base,
      meth_mat = meth_mat,
      extra_mats = extra_mats,
      lo.count = eff_lo_count,
      hi.count = eff_hi_count_internal,
      lo.perc = lo.perc,
      hi.perc = hi.perc,
      mode = mode,
      keep_rule = keep_rule,
      max_missing_prop = max_missing_prop,
      return_what = return_what,
      verbose = FALSE
    )
    
    if (verbose) {
      obj <- out[[i]]
      cov_out <- if (return_what == "all") obj$cov else if (return_what == "cov") obj else NULL
      if (!is.null(cov_out)) {
        n_cpg_nonallna <- sum(rowSums(!is.na(cov_out)) > 0)
        message(
          "  Output summary: ",
          n_cpg_nonallna, " / ", nrow(cov_out),
          " CpGs retained with at least one passing value; ",
          ncol(cov_out), " samples."
        )
      } else {
        n_cpg_nonallna <- sum(rowSums(!is.na(obj)) > 0)
        message(
          "  Output summary: ",
          n_cpg_nonallna, " / ", nrow(obj),
          " CpGs retained with at least one passing value; ",
          ncol(obj), " samples."
        )
      }
    }
  }
  
  out
}
# 


# 
# Align array probes to user-supplied CpG mappings, averaging available probes.
align_array_beta_to_cpgs <- function(beta_matrix, mapping, probe_col = "probe_ids",
                                     cpg_id_col = "cpg_id", sample_metadata = NULL,
                                     array_sample_col = "array_sample_id",
                                     sample_id_col = "sample_id", na.rm = TRUE) {
  beta_matrix <- .to_numeric_matrix(beta_matrix, "beta_matrix")
  if (is.null(rownames(beta_matrix)) || anyDuplicated(rownames(beta_matrix)))
    stop("beta_matrix must have unique probe rownames")
  if (!all(c(probe_col, cpg_id_col) %in% names(mapping)))
    stop("mapping must contain the requested probe and CpG ID columns")
  ids <- as.character(mapping[[cpg_id_col]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Mapping CpG IDs must be unique and nonmissing")
  probes <- strsplit(as.character(mapping[[probe_col]]), ",", fixed=TRUE)
  out <- matrix(NA_real_, nrow=nrow(mapping), ncol=ncol(beta_matrix),
                dimnames=list(ids, colnames(beta_matrix)))
  for (i in seq_along(probes)) {
    idx <- match(unique(trimws(probes[[i]])), rownames(beta_matrix))
    idx <- idx[!is.na(idx)]
    if (length(idx)) out[i, ] <- colMeans(beta_matrix[idx, , drop=FALSE], na.rm=na.rm)
  }
  out[is.nan(out)] <- NA_real_
  if (!is.null(sample_metadata)) {
    if (!all(c(array_sample_col, sample_id_col) %in% names(sample_metadata)))
      stop("Missing sample metadata columns")
    if (anyDuplicated(sample_metadata[[array_sample_col]])) stop("Array sample IDs must be unique")
    idx <- match(colnames(out), sample_metadata[[array_sample_col]])
    names_new <- as.character(sample_metadata[[sample_id_col]][idx])
    if (anyNA(names_new) || any(!nzchar(names_new)) || anyDuplicated(names_new))
      stop("Every matrix sample must map to a unique nonmissing sample ID")
    colnames(out) <- names_new
  }
  out
}

beta_to_m <- function(beta_mat, eps = 1e-6) {
  # Check input
  if (!is.matrix(beta_mat) && !is.data.frame(beta_mat)) {
    stop("beta_mat must be a matrix or data.frame")
  }

  beta_mat <- as.matrix(beta_mat)

  if (!is.numeric(beta_mat)) {
    stop("beta_mat must contain numeric values")
  }

  # Warn if values are outside expected beta range
  if (any(beta_mat < 0 | beta_mat > 1, na.rm = TRUE)) {
    stop("beta values must be between 0 and 1")
  }

  # Clamp values to avoid Inf from log2(0) or division by zero
  beta_mat <- pmax(pmin(beta_mat, 1 - eps), eps)

  # Convert beta values to M values
  m_mat <- log2(beta_mat / (1 - beta_mat))

  return(m_mat)
}




# Takes full beta matrix, returns vector of means
calc_mean_cpg_beta <- function(
    beta_mat,
    na.rm = TRUE,
    min_non_missing = 1
) {
  if (is.data.frame(beta_mat)) beta_mat <- as.matrix(beta_mat)
  if (!is.matrix(beta_mat)) stop("beta_mat must be a matrix or data.frame.")
  if (is.null(rownames(beta_mat))) stop("beta_mat must have CpG/probe rownames.")
  
  suppressWarnings(mode(beta_mat) <- "numeric")
  
  n_non_missing <- rowSums(is.finite(beta_mat))
  
  mean_beta <- rowMeans(beta_mat, na.rm = na.rm)
  
  mean_beta[n_non_missing < min_non_missing] <- NA_real_
  
  out <- mean_beta
  names(out) <- rownames(beta_mat)
  
  out
}


# Classify mean beta values into hypo / hemi / hyper
classify_beta_state <- function(
    beta_vec,
    hypo_cutoff = 0.30,
    hyper_cutoff = 0.70,
    ordered_factor = TRUE
) {
  if (is.null(names(beta_vec))) {
    stop("beta_vec should be a named vector so CpGs can be tracked.")
  }
  
  beta_state <- rep(NA_character_, length(beta_vec))
  
  beta_state[is.finite(beta_vec) & beta_vec < hypo_cutoff] <- "hypo"
  beta_state[is.finite(beta_vec) & beta_vec >= hypo_cutoff & beta_vec <= hyper_cutoff] <- "hemi"
  beta_state[is.finite(beta_vec) & beta_vec > hyper_cutoff] <- "hyper"
  
  names(beta_state) <- names(beta_vec)
  
  if (ordered_factor) {
    beta_state <- factor(
      beta_state,
      levels = c("hypo", "hemi", "hyper"),
      ordered = TRUE
    )
  }
  
  beta_state
}


# Summarize one beta matrix
summarize_cpg_beta_states <- function(
    beta_mat,
    platform_name = "platform",
    hypo_cutoff = 0.30,
    hyper_cutoff = 0.70,
    min_non_missing = 1
) {
  mean_beta <- calc_mean_cpg_beta(
    beta_mat = beta_mat,
    na.rm = TRUE,
    min_non_missing = min_non_missing
  )
  
  beta_state <- classify_beta_state(
    beta_vec = mean_beta,
    hypo_cutoff = hypo_cutoff,
    hyper_cutoff = hyper_cutoff
  )
  
  data.frame(
    cpg_id = names(mean_beta),
    platform = platform_name,
    mean_beta = as.numeric(mean_beta),
    beta_state = beta_state,
    stringsAsFactors = FALSE
  )
}
