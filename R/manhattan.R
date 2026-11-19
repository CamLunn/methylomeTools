plot_dss_manhattan <- function(
    dss_results,
    chr_col = "chr",
    pos_col = "pos",
    p_col = "pvals",
    stat_col = "stat",
    fdr_col = "fdr_global",
    directional = FALSE,
    fdr_threshold = 0.05,
    nominal_p_threshold = NULL,
    point_size = 0.60,
    point_alpha = 0.60,
    chromosome_colors = c("#404040", "#BDBDBD"),
    hypo_color = "#2166AC",
    hyper_color = "#B2182B",
    highlight_fdr = TRUE,
    label_top_n = 0,
    label_col = NULL,
    label_fdr_only = TRUE,
    max_logp = NULL,
    symmetric_y_max = NULL,
    chromosome_gap = 1e6,
    show_chr_boundaries = TRUE,
    chr_boundary_color = "#D9D9D9",
    chr_boundary_linewidth = 0.30,
    include_chr = NULL,
    exclude_chr = NULL,
    title = NULL,
    subtitle = NULL
) {
  
  # -------------------------------------------------------------------------
  # Package checks
  # -------------------------------------------------------------------------
  
  required_packages <- c(
    "dplyr",
    "ggplot2",
    "rlang",
    "scales"
  )
  
  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(missing_packages) > 0) {
    stop(
      "Install the following packages first: ",
      paste(missing_packages, collapse = ", ")
    )
  }
  
  if (label_top_n > 0 &&
      !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Install 'ggrepel' to label CpG sites.")
  }
  
  # -------------------------------------------------------------------------
  # Validate arguments and columns
  # -------------------------------------------------------------------------
  
  required_cols <- c(chr_col, pos_col, p_col)
  
  if (directional) {
    required_cols <- c(required_cols, stat_col)
  }
  
  if (highlight_fdr || label_fdr_only) {
    required_cols <- c(required_cols, fdr_col)
  }
  
  required_cols <- unique(required_cols)
  
  missing_cols <- setdiff(required_cols, names(dss_results))
  
  if (length(missing_cols) > 0) {
    stop(
      "The following columns are missing from dss_results: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.numeric(fdr_threshold) ||
      length(fdr_threshold) != 1 ||
      !is.finite(fdr_threshold) ||
      fdr_threshold <= 0 ||
      fdr_threshold > 1) {
    stop("fdr_threshold must be one numeric value between 0 and 1.")
  }
  
  if (!is.numeric(chromosome_gap) ||
      length(chromosome_gap) != 1 ||
      !is.finite(chromosome_gap) ||
      chromosome_gap < 0) {
    stop("chromosome_gap must be one non-negative numeric value.")
  }
  
  if (length(chromosome_colors) < 2) {
    stop("chromosome_colors must contain at least two colors.")
  }
  
  if (!is.numeric(label_top_n) ||
      length(label_top_n) != 1 ||
      !is.finite(label_top_n) ||
      label_top_n < 0) {
    stop("label_top_n must be a non-negative number.")
  }
  
  label_top_n <- as.integer(label_top_n)
  
  # -------------------------------------------------------------------------
  # Prepare input data
  # -------------------------------------------------------------------------
  
  dat <- dss_results |>
    dplyr::mutate(
      .chr_raw = as.character(.data[[chr_col]]),
      .chr_clean = sub(
        "^chr",
        "",
        .chr_raw,
        ignore.case = TRUE
      ),
      .chr_clean = toupper(.chr_clean),
      .pos = suppressWarnings(
        as.numeric(.data[[pos_col]])
      ),
      .p = suppressWarnings(
        as.numeric(.data[[p_col]])
      ),
      .stat = if (stat_col %in% names(dss_results)) {
        suppressWarnings(
          as.numeric(.data[[stat_col]])
        )
      } else {
        NA_real_
      },
      .fdr = if (fdr_col %in% names(dss_results)) {
        suppressWarnings(
          as.numeric(.data[[fdr_col]])
        )
      } else {
        NA_real_
      }
    )
  
  # Standardize chromosome filters.
  if (!is.null(include_chr)) {
    include_chr <- toupper(
      sub(
        "^chr",
        "",
        as.character(include_chr),
        ignore.case = TRUE
      )
    )
    
    dat <- dat |>
      dplyr::filter(.chr_clean %in% include_chr)
  }
  
  if (!is.null(exclude_chr)) {
    exclude_chr <- toupper(
      sub(
        "^chr",
        "",
        as.character(exclude_chr),
        ignore.case = TRUE
      )
    )
    
    dat <- dat |>
      dplyr::filter(!.chr_clean %in% exclude_chr)
  }
  
  n_input <- nrow(dat)
  
  # Remove rows without usable chromosome, position, or P-value data.
  dat <- dat |>
    dplyr::filter(
      !is.na(.chr_clean),
      .chr_clean != "",
      is.finite(.pos),
      .pos >= 0,
      is.finite(.p),
      .p >= 0,
      .p <= 1
    )
  
  # A statistic is required only for the directional plot.
  if (directional) {
    dat <- dat |>
      dplyr::filter(is.finite(.stat))
  }
  
  n_valid <- nrow(dat)
  
  if (n_valid == 0) {
    stop(
      "No rows contain valid chromosome, position, P-value",
      if (directional) ", and statistic" else "",
      " data."
    )
  }
  
  # -------------------------------------------------------------------------
  # Handle P-values and calculate -log10(P)
  # -------------------------------------------------------------------------
  
  positive_p <- dat$.p[
    is.finite(dat$.p) & dat$.p > 0
  ]
  
  if (length(positive_p) > 0) {
    smallest_positive_p <- min(positive_p)
  } else {
    smallest_positive_p <- .Machine$double.xmin
  }
  
  replacement_p <- max(
    smallest_positive_p / 10,
    .Machine$double.xmin
  )
  
  dat <- dat |>
    dplyr::mutate(
      .p_plot = dplyr::if_else(
        .p == 0,
        replacement_p,
        .p
      ),
      .logp_raw = -log10(.p_plot)
    )
  
  # Optionally cap extremely large -log10(P) values.
  if (is.null(max_logp)) {
    dat <- dat |>
      dplyr::mutate(.logp = .logp_raw)
  } else {
    
    if (!is.numeric(max_logp) ||
        length(max_logp) != 1 ||
        !is.finite(max_logp) ||
        max_logp <= 0) {
      stop("max_logp must be NULL or one positive numeric value.")
    }
    
    dat <- dat |>
      dplyr::mutate(
        .logp = pmin(.logp_raw, max_logp)
      )
  }
  
  # -------------------------------------------------------------------------
  # Natural chromosome ordering
  # -------------------------------------------------------------------------
  
  observed_chr <- unique(dat$.chr_clean)
  
  chr_numeric <- suppressWarnings(
    as.integer(observed_chr)
  )
  
  chr_order_table <- data.frame(
    .chr_clean = observed_chr,
    .chr_group_order = dplyr::case_when(
      !is.na(chr_numeric) ~ 1L,
      observed_chr == "X" ~ 2L,
      observed_chr == "Y" ~ 3L,
      observed_chr %in% c("M", "MT") ~ 4L,
      TRUE ~ 5L
    ),
    .chr_number = dplyr::case_when(
      !is.na(chr_numeric) ~ chr_numeric,
      observed_chr == "X" ~ 23L,
      observed_chr == "Y" ~ 24L,
      observed_chr %in% c("M", "MT") ~ 25L,
      TRUE ~ 999L
    ),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(
      .chr_group_order,
      .chr_number,
      .chr_clean
    ) |>
    dplyr::mutate(
      .chr_index = dplyr::row_number()
    )
  
  dat <- dat |>
    dplyr::inner_join(
      chr_order_table,
      by = ".chr_clean"
    ) |>
    dplyr::arrange(.chr_index, .pos)
  
  # -------------------------------------------------------------------------
  # Calculate cumulative genomic positions
  # -------------------------------------------------------------------------
  
  chr_info <- dat |>
    dplyr::group_by(.chr_clean, .chr_index) |>
    dplyr::summarise(
      .chr_min = min(.pos),
      .chr_max = max(.pos),
      .chr_length = .chr_max - .chr_min + 1,
      .groups = "drop"
    ) |>
    dplyr::arrange(.chr_index) |>
    dplyr::mutate(
      .offset = dplyr::lag(
        cumsum(.chr_length + chromosome_gap),
        default = 0
      ),
      .plot_start = .offset,
      .plot_end = .offset + .chr_length,
      .center = .offset + (.chr_length / 2)
    )
  
  dat <- dat |>
    dplyr::left_join(
      dplyr::select(
        chr_info,
        .chr_clean,
        .chr_index,
        .chr_min,
        .chr_max,
        .chr_length,
        .offset,
        .plot_start,
        .plot_end,
        .center
      ),
      by = c(".chr_clean", ".chr_index")
    ) |>
    dplyr::mutate(
      .cumulative_pos = .offset + (.pos - .chr_min),
      .chromosome_group = factor(.chr_index %% 2),
      .significant = (
        !is.na(.fdr) &
          is.finite(.fdr) &
          .fdr < fdr_threshold
      ),
      .direction = dplyr::case_when(
        .stat > 0 ~ "Hyper",
        .stat < 0 ~ "Hypo",
        TRUE ~ "No direction"
      ),
      .plot_y = if (directional) {
        sign(.stat) * .logp
      } else {
        .logp
      }
    )
  
  # Boundaries are placed halfway across the gap between chromosomes.
  if (nrow(chr_info) > 1) {
    chromosome_boundaries <- (
      chr_info$.plot_end[-nrow(chr_info)] +
        chr_info$.plot_start[-1]
    ) / 2
  } else {
    chromosome_boundaries <- numeric(0)
  }
  
  # -------------------------------------------------------------------------
  # Symmetric directional y-axis
  # -------------------------------------------------------------------------
  
  if (directional) {
    
    if (is.null(symmetric_y_max)) {
      y_axis_max <- max(
        abs(dat$.plot_y),
        na.rm = TRUE
      )
      
      # Round upward to a readable interval of five.
      y_axis_max <- ceiling(y_axis_max / 5) * 5
    } else {
      
      if (!is.numeric(symmetric_y_max) ||
          length(symmetric_y_max) != 1 ||
          !is.finite(symmetric_y_max) ||
          symmetric_y_max <= 0) {
        stop(
          "symmetric_y_max must be NULL or one positive ",
          "numeric value."
        )
      }
      
      y_axis_max <- symmetric_y_max
    }
    
    if (!is.finite(y_axis_max) || y_axis_max <= 0) {
      y_axis_max <- 5
    }
    
    positive_breaks <- pretty(c(0, y_axis_max))
    
    positive_breaks <- positive_breaks[
      positive_breaks >= 0 &
        positive_breaks <= y_axis_max
    ]
    
    y_breaks <- sort(
      unique(c(-positive_breaks, positive_breaks))
    )
  }
  
  # -------------------------------------------------------------------------
  # X-axis labels
  # -------------------------------------------------------------------------
  
  axis_labels <- chr_info$.chr_clean
  axis_labels[axis_labels == "MT"] <- "M"
  
  
  # -------------------------------------------------------------------------
  # Calculate coutns of significant CpGs
  # -------------------------------------------------------------------------
  n_hyper_fdr <- sum(
    dat$.significant & dat$.direction == "Hyper",
    na.rm = TRUE
  )
  
  n_hypo_fdr <- sum(
    dat$.significant & dat$.direction == "Hypo",
    na.rm = TRUE
  )
  
  # -------------------------------------------------------------------------
  # Construct base plot
  # -------------------------------------------------------------------------
  
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .cumulative_pos,
      y = .plot_y
    )
  )
  
  # Light chromosome separators are drawn behind the data points.
  if (show_chr_boundaries &&
      length(chromosome_boundaries) > 0) {
    
    p <- p +
      ggplot2::geom_vline(
        xintercept = chromosome_boundaries,
        color = chr_boundary_color,
        linewidth = chr_boundary_linewidth,
        linetype = "solid"
      )
  }
  
  # Alternating chromosome colors use fill, leaving color available
  # for significant hyper/hypomethylated loci.
  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(color = .chromosome_group),
      shape = 16,
      size = point_size,
      alpha = point_alpha,
      stroke = 0
    ) +
    ggplot2::scale_x_continuous(
      breaks = chr_info$.center,
      labels = axis_labels,
      expand = ggplot2::expansion(
        mult = c(0.005, 0.005)
      )
    ) +
    ggplot2::labs(
      x = "Chromosome",
      y = if (directional) {
        expression("Signed " * -log[10](italic(P)))
      } else {
        expression(-log[10](italic(P)))
      },
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(
        color = "black",
        linewidth = 0.45
      ),
      axis.line.y = ggplot2::element_line(
        color = "black",
        linewidth = 0.45
      ),
      axis.ticks.x = ggplot2::element_line(
        color = "grey35",
        linewidth = 0.35
      ),
      axis.ticks.y = ggplot2::element_line(
        color = "black",
        linewidth = 0.35
      ),
      axis.text.x = ggplot2::element_text(
        color = "black",
        size = 9,
        margin = ggplot2::margin(t = 4)
      ),
      axis.text.y = ggplot2::element_text(
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        color = "black"
      ),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 7)
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 7)
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      plot.subtitle = ggplot2::element_text(
        color = "grey30",
        size = 10
      ),
      legend.position = "top",
      legend.justification = "left",
      legend.text = ggplot2::element_text(
        size = 9
      ),
      plot.margin = ggplot2::margin(
        t = 8,
        r = 12,
        b = 8,
        l = 8
      )
    )
  
  # -------------------------------------------------------------------------
  # Y-axis scales
  # -------------------------------------------------------------------------
  
  if (directional) {
    p <- p +
      ggplot2::scale_y_continuous(
        limits = c(-y_axis_max, y_axis_max),
        breaks = y_breaks,
        labels = function(x) {
          scales::label_number(accuracy = 1)(abs(x))
        },
        expand = ggplot2::expansion(
          mult = c(0.02, 0.02)
        )
      )
  } else {
    p <- p +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(
          mult = c(0, 0.06)
        ),
        labels = scales::label_number(
          accuracy = 1
        )
      )
  }
  
  # Central x-axis in directional mode.
  if (directional) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = 0,
        color = "black",
        linewidth = 0.50
      )
  }
  
  # -------------------------------------------------------------------------
  # Highlight FDR-significant loci
  # -------------------------------------------------------------------------
  
  if (highlight_fdr &&
      any(dat$.significant, na.rm = TRUE)) {
    
    if (directional) {
      
      significant_dat <- dat |>
        dplyr::filter(
          .significant,
          .direction %in% c("Hyper", "Hypo")
        )
      
      if (nrow(significant_dat) > 0) {
        p <- p +
          ggplot2::geom_point(
            data = significant_dat,
            ggplot2::aes(color = .direction),
            shape = 16,
            size = max(point_size * 1.5, 0.65),
            alpha = 0.95,
            stroke = 0
          )
      }
      
    } else {
      
      p <- p +
        ggplot2::geom_point(
          data = dplyr::filter(dat, .significant),
          color = hyper_color,
          shape = 16,
          size = max(point_size * 1.5, 0.65),
          alpha = 0.95,
          stroke = 0
        )
    }
  }
  
  
  if (directional && highlight_fdr) {
    
    p <- p +
      ggplot2::scale_color_manual(
        name = NULL,
        values = c(
          "0" = chromosome_colors[1],
          "1" = chromosome_colors[2],
          "Hyper" = hyper_color,
          "Hypo" = hypo_color
        ),
        breaks = c("Hyper", "Hypo"),
        labels = c(
          paste0(
            "Hypermethylated (n = ",
            scales::comma(n_hyper_fdr),
            ")"
          ),
          paste0(
            "Hypomethylated (n = ",
            scales::comma(n_hypo_fdr),
            ")"
          )
        ),
        guide = ggplot2::guide_legend(
          override.aes = list(
            size = 3,
            alpha = 1,
            shape = 16
          )
        )
      )
    
  } else {
    
    p <- p +
      ggplot2::scale_color_manual(
        values = c(
          "0" = chromosome_colors[1],
          "1" = chromosome_colors[2]
        ),
        guide = "none"
      )
  }
  
  # -------------------------------------------------------------------------
  # Nominal raw-P threshold
  # -------------------------------------------------------------------------
  
  if (!is.null(nominal_p_threshold)) {
    
    if (!is.numeric(nominal_p_threshold) ||
        length(nominal_p_threshold) != 1 ||
        !is.finite(nominal_p_threshold) ||
        nominal_p_threshold <= 0 ||
        nominal_p_threshold > 1) {
      stop(
        "nominal_p_threshold must be one numeric value ",
        "between 0 and 1."
      )
    }
    
    threshold_y <- -log10(nominal_p_threshold)
    
    p <- p +
      ggplot2::geom_hline(
        yintercept = if (directional) {
          c(-threshold_y, threshold_y)
        } else {
          threshold_y
        },
        linetype = "dashed",
        linewidth = 0.45,
        color = "#525252"
      )
  }
  
  # -------------------------------------------------------------------------
  # Label strongest loci
  # -------------------------------------------------------------------------
  
  label_data <- NULL
  
  if (label_top_n > 0) {
    
    label_candidates <- dat
    
    if (label_fdr_only) {
      label_candidates <- label_candidates |>
        dplyr::filter(.significant)
    }
    
    if (nrow(label_candidates) > 0) {
      
      label_data <- label_candidates |>
        dplyr::arrange(.p) |>
        dplyr::slice_head(n = label_top_n)
      
      if (is.null(label_col)) {
        
        label_data <- label_data |>
          dplyr::mutate(
            .plot_label = paste0(
              "chr",
              .chr_clean,
              ":",
              format(
                .pos,
                scientific = FALSE,
                trim = TRUE
              )
            )
          )
        
      } else {
        
        if (!label_col %in% names(label_data)) {
          stop(
            "label_col '",
            label_col,
            "' is not present in dss_results."
          )
        }
        
        label_data$.plot_label <- as.character(
          label_data[[label_col]]
        )
      }
      
      p <- p +
        ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(label = .plot_label),
          size = 3,
          color = "black",
          box.padding = 0.35,
          point.padding = 0.20,
          min.segment.length = 0,
          segment.color = "grey50",
          segment.linewidth = 0.30,
          max.overlaps = Inf,
          seed = 1,
          show.legend = FALSE
        )
    }
  }
  
  # -------------------------------------------------------------------------
  # Diagnostics
  # -------------------------------------------------------------------------
  
  diagnostics <- data.frame(
    n_input = n_input,
    n_plotted = n_valid,
    n_excluded = n_input - n_valid,
    n_zero_p_replaced = sum(dat$.p == 0, na.rm = TRUE),
    n_fdr_significant = sum(
      dat$.significant,
      na.rm = TRUE
    ),
    n_hyper_fdr = sum(
      dat$.significant &
        dat$.direction == "Hyper",
      na.rm = TRUE
    ),
    n_hypo_fdr = sum(
      dat$.significant &
        dat$.direction == "Hypo",
      na.rm = TRUE
    ),
    replacement_p = replacement_p,
    directional = directional,
    symmetric_y_max = if (directional) {
      y_axis_max
    } else {
      NA_real_
    }
  )
  
  structure(
    list(
      plot = p,
      plot_data = dat,
      chromosome_info = chr_info,
      chromosome_boundaries = chromosome_boundaries,
      labeled_loci = label_data,
      diagnostics = diagnostics
    ),
    class = "dss_manhattan_result"
  )
}
