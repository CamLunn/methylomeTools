.extract_top_dml_betas <- function(
    bs,
    dml_table,
    meta,
    sample_id_col = "sample_id",
    group_col,
    fdr_col = "fdr_global",
    chr_col = "chr",
    pos_col = "pos",
    top_n = 9,
    cpg_ids = NULL,
    min_cov = 1
) {
  stopifnot(inherits(bs, "BSseq"))
  stopifnot(is.data.frame(dml_table))
  stopifnot(is.data.frame(meta))
  
  required_dml_cols <- c(chr_col, pos_col, fdr_col)
  missing_dml_cols <- setdiff(required_dml_cols, colnames(dml_table))
  
  if (length(missing_dml_cols) > 0) {
    stop(
      "Missing columns in dml_table: ",
      paste(missing_dml_cols, collapse = ", ")
    )
  }
  
  required_meta_cols <- c(sample_id_col, group_col)
  missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
  
  if (length(missing_meta_cols) > 0) {
    stop(
      "Missing columns in meta: ",
      paste(missing_meta_cols, collapse = ", ")
    )
  }
  
  assay_names <- SummarizedExperiment::assayNames(bs)
  
  if (!all(c("M", "Cov") %in% assay_names)) {
    stop(
      "Expected BSseq assays named 'M' and 'Cov'. Found assays: ",
      paste(assay_names, collapse = ", ")
    )
  }
  
  # Accept either:
  #   "chr1.12345"
  # or
  #   "chr1.12345.12345"
  #
  # Internally, all IDs are standardized as "chr.pos".
  normalize_cpg_ids <- function(ids) {
    ids <- trimws(as.character(ids))
    
    sub(
      pattern = "^(.+)\\.([0-9]+)(?:\\.[0-9]+)?$",
      replacement = "\\1.\\2",
      x = ids
    )
  }
  
  # Standardize DML table with CpG IDs in chr.pos format
  dml_table2 <- dml_table %>%
    dplyr::mutate(
      .dml_chr = as.character(.data[[chr_col]]),
      .dml_pos = as.integer(.data[[pos_col]]),
      cpg_id = paste0(.dml_chr, ".", .dml_pos),
      cpg_label = paste0(
        .dml_chr,
        ":",
        .dml_pos,
        "\nFDR=",
        signif(.data[[fdr_col]], 3)
      )
    )
  
  # ----------------------------------------------------------
  # Select CpGs:
  # - If cpg_ids is supplied, use those requested CpGs.
  # - Otherwise, select the top_n CpGs by FDR.
  # ----------------------------------------------------------
  if (!is.null(cpg_ids)) {
    cpg_ids <- unique(normalize_cpg_ids(cpg_ids))
    
    if (length(cpg_ids) == 0) {
      stop("cpg_ids was provided but is empty.")
    }
    
    top_dml <- dml_table2 %>%
      dplyr::filter(cpg_id %in% cpg_ids) %>%
      dplyr::mutate(.cpg_order = match(cpg_id, cpg_ids)) %>%
      dplyr::arrange(.cpg_order) %>%
      dplyr::select(-.cpg_order)
    
    missing_from_dml <- setdiff(cpg_ids, top_dml$cpg_id)
    
    if (length(missing_from_dml) > 0) {
      warning(
        length(missing_from_dml),
        " requested cpg_ids were not found in dml_table and will be dropped: ",
        paste(missing_from_dml, collapse = ", ")
      )
    }
    
    if (nrow(top_dml) == 0) {
      stop("None of the requested cpg_ids were found in dml_table.")
    }
    
  } else {
    top_dml <- dml_table2 %>%
      dplyr::filter(!is.na(.data[[fdr_col]])) %>%
      dplyr::arrange(.data[[fdr_col]]) %>%
      dplyr::slice_head(n = top_n)
    
    if (nrow(top_dml) == 0) {
      stop(
        "No valid DMLs found after filtering non-missing ",
        fdr_col,
        "."
      )
    }
  }
  
  # ----------------------------------------------------------
  # Match DML positions to rows in the BSseq object
  # ----------------------------------------------------------
  gr <- GenomicRanges::granges(bs)
  
  bs_loci <- tibble::tibble(
    bs_row = seq_along(gr),
    .dml_chr = as.character(GenomeInfoDb::seqnames(gr)),
    .dml_pos = as.integer(GenomicRanges::start(gr))
  )
  
  top_dml_matched <- top_dml %>%
    dplyr::left_join(
      bs_loci,
      by = c(".dml_chr", ".dml_pos")
    )
  
  missing_loci <- top_dml_matched %>%
    dplyr::filter(is.na(bs_row))
  
  if (nrow(missing_loci) > 0) {
    warning(
      nrow(missing_loci),
      " selected CpGs were not found in the BSseq object and will be dropped: ",
      paste(missing_loci$cpg_id, collapse = ", ")
    )
  }
  
  top_dml_matched <- top_dml_matched %>%
    dplyr::filter(!is.na(bs_row))
  
  if (nrow(top_dml_matched) == 0) {
    stop("None of the selected CpGs were found in the BSseq object.")
  }
  
  # Preserve explicitly requested CpG order.
  # For top_n, preserve FDR-ranked order.
  if (!is.null(cpg_ids)) {
    top_dml_matched <- top_dml_matched %>%
      dplyr::mutate(.cpg_order = match(cpg_id, cpg_ids)) %>%
      dplyr::arrange(.cpg_order) %>%
      dplyr::select(-.cpg_order)
  }
  
  # ----------------------------------------------------------
  # Extract methylated-read count (M) and coverage (Cov)
  # ----------------------------------------------------------
  meth_mat <- SummarizedExperiment::assay(
    bs,
    "M"
  )[top_dml_matched$bs_row, , drop = FALSE]
  
  cov_mat <- SummarizedExperiment::assay(
    bs,
    "Cov"
  )[top_dml_matched$bs_row, , drop = FALSE]
  
  meth_mat <- as.matrix(meth_mat)
  cov_mat <- as.matrix(cov_mat)
  
  storage.mode(meth_mat) <- "numeric"
  storage.mode(cov_mat) <- "numeric"
  
  beta_mat <- meth_mat / cov_mat
  
  # Mark sites below the minimum coverage threshold as missing.
  beta_mat[cov_mat < min_cov] <- NA
  beta_mat[!is.finite(beta_mat)] <- NA
  
  rownames(beta_mat) <- top_dml_matched$cpg_id
  rownames(cov_mat) <- top_dml_matched$cpg_id
  rownames(meth_mat) <- top_dml_matched$cpg_id
  
  # ----------------------------------------------------------
  # Match BSseq samples to metadata
  # ----------------------------------------------------------
  sample_ids_bs <- colnames(beta_mat)
  
  meta_clean <- meta %>%
    dplyr::mutate(
      .sample_id = as.character(.data[[sample_id_col]]),
      .group_value = .data[[group_col]]
    ) %>%
    dplyr::filter(.sample_id %in% sample_ids_bs)
  
  if (nrow(meta_clean) == 0) {
    stop(
      "No samples in meta matched colnames(bs). ",
      "Check sample_id_col and BSseq sample names."
    )
  }
  
  beta_mat <- beta_mat[, meta_clean$.sample_id, drop = FALSE]
  cov_mat <- cov_mat[, meta_clean$.sample_id, drop = FALSE]
  meth_mat <- meth_mat[, meta_clean$.sample_id, drop = FALSE]
  
  # ----------------------------------------------------------
  # Convert matrices to long plotting format
  # ----------------------------------------------------------
  beta_long <- as.data.frame(beta_mat) %>%
    tibble::rownames_to_column("cpg_id") %>%
    tidyr::pivot_longer(
      cols = -cpg_id,
      names_to = ".sample_id",
      values_to = "beta"
    ) %>%
    dplyr::left_join(
      as.data.frame(cov_mat) %>%
        tibble::rownames_to_column("cpg_id") %>%
        tidyr::pivot_longer(
          cols = -cpg_id,
          names_to = ".sample_id",
          values_to = "coverage"
        ),
      by = c("cpg_id", ".sample_id")
    ) %>%
    dplyr::left_join(
      as.data.frame(meth_mat) %>%
        tibble::rownames_to_column("cpg_id") %>%
        tidyr::pivot_longer(
          cols = -cpg_id,
          names_to = ".sample_id",
          values_to = "meth_count"
        ),
      by = c("cpg_id", ".sample_id")
    ) %>%
    dplyr::left_join(
      meta_clean,
      by = ".sample_id"
    ) %>%
    dplyr::left_join(
      top_dml_matched %>%
        dplyr::select(
          cpg_id,
          cpg_label,
          dplyr::all_of(fdr_col),
          dplyr::all_of(chr_col),
          dplyr::all_of(pos_col)
        ),
      by = "cpg_id"
    )
  
  # Keep facet order stable and aligned to the selected/ranked CpGs.
  facet_levels <- top_dml_matched$cpg_label
  plot_levels <- top_dml_matched$cpg_id
  
  beta_long <- beta_long %>%
    dplyr::mutate(
      cpg_id = factor(cpg_id, levels = plot_levels),
      cpg_label = factor(cpg_label, levels = facet_levels)
    )
  
  list(
    beta_long = beta_long,
    top_dml = top_dml_matched,
    beta_mat = beta_mat,
    cov_mat = cov_mat,
    meth_mat = meth_mat
  )
}


# ------------------------------------------------------------
# Function 1: Boxplots using tertiles of a continuous variable
# ------------------------------------------------------------
plot_top_dml_betas_boxplot <- function(
    bs,
    dml_table,
    meta,
    sample_id_col = "sample_id",
    group_col,
    fdr_col = "fdr_global",
    chr_col = "chr",
    pos_col = "pos",
    top_n = 9,
    cpg_ids = NULL,
    min_cov = 1,
    tertile_labels = c("Low", "Middle", "High"),
    pdf_file = NULL,
    ncol = 3,
    point_alpha = 0.75,
    point_size = 1.8
) {
  extracted <- .extract_top_dml_betas(
    bs = bs,
    dml_table = dml_table,
    meta = meta,
    sample_id_col = sample_id_col,
    group_col = group_col,
    fdr_col = fdr_col,
    chr_col = chr_col,
    pos_col = pos_col,
    top_n = top_n,
    cpg_ids = cpg_ids,
    min_cov = min_cov
  )
  
  plot_df <- extracted$beta_long %>%
    dplyr::filter(!is.na(.group_value), !is.na(beta)) %>%
    dplyr::mutate(
      group_tertile_num = dplyr::ntile(.group_value, 3),
      group_tertile = factor(
        group_tertile_num,
        levels = 1:3,
        labels = tertile_labels
      )
    )
  
  title_text <- if (!is.null(cpg_ids)) {
    paste0(
      "Selected DML CpGs (n=",
      length(unique(as.character(cpg_ids))),
      ")"
    )
  } else {
    paste0("Top ", top_n, " DML CpGs by ", fdr_col)
  }
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = group_tertile, y = beta)
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      width = 0.65
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      height = 0,
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::facet_wrap(
      ~cpg_label,
      scales = "free_y",
      ncol = ncol
    ) +
    ggplot2::labs(
      x = paste0(group_col, " tertile"),
      y = "Beta value",
      title = title_text,
      subtitle = "Continuous variable grouped into tertiles"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  if (!is.null(pdf_file)) {
    dir.create(
      dirname(pdf_file),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    ggplot2::ggsave(
      filename = pdf_file,
      plot = p,
      width = 11,
      height = 8.5,
      units = "in"
    )
  }
  
  list(
    plot = p,
    plot_df = plot_df,
    extracted = extracted
  )
}


# ------------------------------------------------------------
# Function 2: Scatterplots using the continuous variable directly
# ------------------------------------------------------------
plot_top_dml_betas_scatter <- function(
    bs,
    dml_table,
    meta,
    sample_id_col = "sample_id",
    group_col,
    color_col = NULL,
    fdr_col = "fdr_global",
    chr_col = "chr",
    pos_col = "pos",
    top_n = 9,
    cpg_ids = NULL,
    min_cov = 1,
    pdf_file = NULL,
    ncol = 3,
    point_alpha = 0.75,
    point_size = 1.8,
    add_lm = TRUE
) {
  extracted <- .extract_top_dml_betas(
    bs = bs,
    dml_table = dml_table,
    meta = meta,
    sample_id_col = sample_id_col,
    group_col = group_col,
    fdr_col = fdr_col,
    chr_col = chr_col,
    pos_col = pos_col,
    top_n = top_n,
    cpg_ids = cpg_ids,
    min_cov = min_cov
  )
  
  plot_df <- extracted$beta_long %>%
    dplyr::filter(!is.na(.group_value), !is.na(beta)) %>%
    dplyr::mutate(
      .group_value = as.numeric(.group_value)
    )
  
  if (!is.null(color_col) && !color_col %in% colnames(plot_df)) {
    stop(
      sprintf(
        "color_col '%s' not found in plotting data.",
        color_col
      )
    )
  }
  
  cor_df <- plot_df %>%
    dplyr::group_by(cpg_id, cpg_label) %>%
    dplyr::summarize(
      n = sum(!is.na(beta) & !is.na(.group_value)),
      pearson_r = suppressWarnings(
        stats::cor(
          beta,
          .group_value,
          use = "complete.obs",
          method = "pearson"
        )
      ),
      spearman_r = suppressWarnings(
        stats::cor(
          beta,
          .group_value,
          use = "complete.obs",
          method = "spearman"
        )
      ),
      .groups = "drop"
    )
  
  aes_args <- list(
    x = quote(.group_value),
    y = quote(beta)
  )
  
  if (!is.null(color_col)) {
    aes_args$color <- rlang::sym(color_col)
  }
  
  title_text <- if (!is.null(cpg_ids)) {
    paste0(
      "Selected DML CpGs (n=",
      length(unique(as.character(cpg_ids))),
      ")"
    )
  } else {
    paste0("Top ", top_n, " DML CpGs by ", fdr_col)
  }
  
  p <- ggplot2::ggplot(
    plot_df,
    do.call(ggplot2::aes, aes_args)
  ) +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::facet_wrap(
      ~cpg_label,
      scales = "free_y",
      ncol = ncol
    ) +
    ggplot2::labs(
      x = group_col,
      y = "Beta value",
      color = color_col,
      title = title_text,
      subtitle = "Continuous exposure plotted directly"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 9)
    )
  
  if (isTRUE(add_lm)) {
    p <- p +
      ggplot2::geom_smooth(
        method = "lm",
        se = TRUE,
        linewidth = 0.7
      )
  }
  
  if (!is.null(pdf_file)) {
    dir.create(
      dirname(pdf_file),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    ggplot2::ggsave(
      filename = pdf_file,
      plot = p,
      width = 11,
      height = 8.5,
      units = "in"
    )
  }
  
  list(
    plot = p,
    plot_df = plot_df,
    cor_df = cor_df,
    extracted = extracted
  )
}

