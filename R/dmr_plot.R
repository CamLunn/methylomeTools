plot_dmr_browser <- function(
    dmr_df,
    methylation_long,
    dmr_id,
    phenotype_col,
    
    # Optional annotation RDS paths
    regulatory_df = NULL,
    cpg_islands_df = NULL,
    gene_features_df = NULL,
    
    # Column names
    dmr_id_col = "dmr_id",
    sample_col = "sample_id",
    chr_col = "chr",
    pos_col = "pos",
    beta_col = "beta",
    coverage_col = "coverage",
    
    # Annotation column names
    regulatory_chr_col = "chr",
    regulatory_start_col = "start",
    regulatory_end_col = "end",
    regulatory_type_col = "type",
    
    cpg_chr_col = "chr",
    cpg_start_col = "start",
    cpg_end_col = "end",
    cpg_name_col = NULL,
    
    gene_chr_col = "chr",
    gene_start_col = "start",
    gene_end_col = "end",
    gene_name_col = "gene_name",
    gene_id_col = "gene_id",
    transcript_id_col = "transcript_id",
    gene_feature_col = "feature",
    gene_strand_col = "strand",
    canonical_transcript_id_col = "canonical_transcript_id",
    
    # Window
    region = NULL,
    padding = 500,
    
    # Phenotype display
    phenotype_type = c(
      "auto",
      "categorical",
      "continuous"
    ),
    continuous_display = c(
      "continuous",
      "tertiles"
    ),
    group_colors = NULL,
    continuous_palette = c(
      "#313695",
      "#74ADD1",
      "#FFFFBF",
      "#F46D43",
      "#A50026"
    ),
    
    # Methylation appearance
    plot_means = FALSE,
    connect_means = TRUE,
    mean_linewidth = 0.75,
    curve_method = c(
      "loess",
      "gam",
      "none"
    ),
    curve_span = 0.70,
    curve_alpha = 0.55,
    curve_linewidth = 0.55,
    point_alpha = 0.80,
    point_size_range = c(.5, 3),
    coverage_limits = NULL,
    
    # Annotation appearance
    regulatory_colors = c(
      promoter = "#D73027",
      enhancer = "#FFD92F",
      ctcf = "#80CBE5"
    ),
    cpg_island_color = "#31A354",
    exon_color = "#3B3B3B",
    intron_color = "#696969",
    tss_color = "#7B3294",
    dmr_highlight_color = "#FF8C00",
    dmr_highlight_alpha = 0.13,
    tss_arrow_length_prop = 0.025,
    tss_vertical_nudge = 0.18,
    intron_arrow_spacing = NULL,
    intron_arrow_length_prop = 0.35,
    annotation_arrow_size = 0.12,
    
    # Labels and layout
    title = NULL,
    show_sample_legend = TRUE,
    max_gene_transcripts = 20,
    relative_track_heights = c(
      methylation = 5.5,
      genes = 2.2,
      cpg = 0.65,
      regulatory = 1
    ),
    base_size = 11
) {
  
  # ============================================================
  # Dependencies
  # ============================================================
  
  required_packages <- c(
    "dplyr",
    "ggplot2",
    "scales",
    "patchwork"
  )
  
  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]
  
  if (length(missing_packages) > 0L) {
    stop(
      "Install the following required packages: ",
      paste(missing_packages, collapse = ", ")
    )
  }
  
  phenotype_type <- match.arg(phenotype_type)
  continuous_display <- match.arg(continuous_display)
  curve_method <- match.arg(curve_method)
  
  if (
    length(plot_means) != 1L ||
    !is.logical(plot_means) ||
    is.na(plot_means)
  ) {
    stop("`plot_means` must be TRUE or FALSE.")
  }
  
  if (
    length(connect_means) != 1L ||
    !is.logical(connect_means) ||
    is.na(connect_means)
  ) {
    stop("`connect_means` must be TRUE or FALSE.")
  }
  
  # ============================================================
  # Read annotation RDS files
  # ============================================================
  
  read_annotation_rds <- function(x, argument_name) {
    
    if (is.null(x)) {
      return(NULL)
    }
    
    if (!is.character(x) || length(x) != 1L || is.na(x)) {
      stop(
        "`", argument_name,
        "` must be NULL or a single path to an RDS file."
      )
    }
    
    if (!file.exists(x)) {
      stop(
        "The file supplied to `",
        argument_name,
        "` does not exist:\n",
        x
      )
    }
    
    object <- readRDS(x)
    
    if (!is.data.frame(object)) {
      stop(
        "The RDS object loaded from `",
        argument_name,
        "` is not a data frame."
      )
    }
    
    as.data.frame(object)
  }
  
  regulatory_data <- read_annotation_rds(
    regulatory_df,
    "regulatory_df"
  )
  
  cpg_islands_data <- read_annotation_rds(
    cpg_islands_df,
    "cpg_islands_df"
  )
  
  gene_features_data <- read_annotation_rds(
    gene_features_df,
    "gene_features_df"
  )
  
  # ============================================================
  # Validate primary inputs
  # ============================================================
  
  dmr_df <- as.data.frame(dmr_df)
  methylation_long <- as.data.frame(methylation_long)
  
  required_dmr_columns <- c(
    dmr_id_col,
    chr_col,
    "start",
    "end"
  )
  
  missing_dmr_columns <- setdiff(
    required_dmr_columns,
    names(dmr_df)
  )
  
  if (length(missing_dmr_columns) > 0L) {
    stop(
      "`dmr_df` is missing required columns: ",
      paste(missing_dmr_columns, collapse = ", ")
    )
  }
  
  required_methylation_columns <- c(
    sample_col,
    chr_col,
    pos_col,
    beta_col,
    coverage_col,
    phenotype_col
  )
  
  missing_methylation_columns <- setdiff(
    required_methylation_columns,
    names(methylation_long)
  )
  
  if (length(missing_methylation_columns) > 0L) {
    stop(
      "`methylation_long` is missing required columns: ",
      paste(missing_methylation_columns, collapse = ", ")
    )
  }
  
  target_dmr_id <- as.character(dmr_id)
  
  if (
    length(target_dmr_id) != 1L ||
    is.na(target_dmr_id) ||
    !nzchar(target_dmr_id)
  ) {
    stop("`dmr_id` must be one nonmissing DMR identifier.")
  }
  
  # Explicitly separate the function argument from the data column.
  selected_dmr <- dmr_df |>
    dplyr::filter(
      as.character(.data[[dmr_id_col]]) ==
        .env$target_dmr_id
    ) |>
    dplyr::transmute(
      dmr_id = as.character(.data[[dmr_id_col]]),
      chr = as.character(.data[[chr_col]]),
      start = as.integer(.data$start),
      end = as.integer(.data$end)
    ) |>
    dplyr::distinct()
  
  if (nrow(selected_dmr) == 0L) {
    stop(
      "`dmr_id = \"",
      target_dmr_id,
      "\"` was not found in `dmr_df`."
    )
  }
  
  if (nrow(selected_dmr) != 1L) {
    stop(
      "`dmr_id = \"",
      target_dmr_id,
      "\"` corresponds to ",
      nrow(selected_dmr),
      " distinct genomic intervals."
    )
  }
  
  dmr_chr <- selected_dmr$chr[[1L]]
  dmr_start <- selected_dmr$start[[1L]]
  dmr_end <- selected_dmr$end[[1L]]
  
  if (
    is.na(dmr_start) ||
    is.na(dmr_end) ||
    dmr_start > dmr_end
  ) {
    stop("The selected DMR has invalid start or end coordinates.")
  }
  
  # ============================================================
  # Determine plotting region
  # ============================================================
  
  if (is.null(region)) {
    region_start <- max(1L, dmr_start - as.integer(padding))
    region_end <- dmr_end + as.integer(padding)
  } else {
    
    if (!is.numeric(region) || length(region) != 2L) {
      stop(
        "`region` must be NULL or a numeric vector: ",
        "c(start, end)."
      )
    }
    
    region_start <- as.integer(region[[1L]])
    region_end <- as.integer(region[[2L]])
    
    if (
      is.na(region_start) ||
      is.na(region_end) ||
      region_start >= region_end
    ) {
      stop("`region` must contain a valid start and end.")
    }
  }
  
  # ============================================================
  # Prepare methylation observations
  # ============================================================
  
  plot_data <- methylation_long |>
    dplyr::mutate(
      "{chr_col}" := as.character(.data[[chr_col]]),
      "{pos_col}" := as.integer(.data[[pos_col]]),
      "{beta_col}" := as.numeric(.data[[beta_col]]),
      "{coverage_col}" := as.numeric(.data[[coverage_col]])
    ) |>
    dplyr::filter(
      .data[[chr_col]] == .env$dmr_chr,
      .data[[pos_col]] >= .env$region_start,
      .data[[pos_col]] <= .env$region_end,
      is.finite(.data[[beta_col]]),
      is.finite(.data[[coverage_col]]),
      !is.na(.data[[phenotype_col]])
    )
  
  if (nrow(plot_data) == 0L) {
    stop(
      "No observed methylation measurements were found in ",
      dmr_chr, ":",
      region_start, "-", region_end, "."
    )
  }
  
  # ============================================================
  # Determine phenotype type
  # ============================================================
  
  if (phenotype_type == "auto") {
    phenotype_vector <- plot_data[[phenotype_col]]
    
    if (
      is.numeric(phenotype_vector) &&
      dplyr::n_distinct(
        phenotype_vector,
        na.rm = TRUE
      ) > 6L
    ) {
      phenotype_type <- "continuous"
    } else {
      phenotype_type <- "categorical"
    }
  }
  
  if (plot_means && phenotype_type != "categorical") {
    stop(
      "`plot_means = TRUE` requires a categorical phenotype. ",
      "Use `phenotype_type = 'categorical'`."
    )
  }
  
  # ============================================================
  # Prepare phenotype display
  # ============================================================
  
  phenotype_display_col <- ".phenotype_display"
  
  if (phenotype_type == "continuous") {
    
    plot_data[[phenotype_col]] <- as.numeric(
      plot_data[[phenotype_col]]
    )
    
    plot_data <- plot_data |>
      dplyr::filter(
        is.finite(.data[[phenotype_col]])
      )
    
    if (continuous_display == "tertiles") {
      
      phenotype_values <- plot_data[[phenotype_col]]
      
      breaks <- unique(
        stats::quantile(
          phenotype_values,
          probs = c(0, 1 / 3, 2 / 3, 1),
          na.rm = TRUE,
          names = FALSE
        )
      )
      
      if (length(breaks) < 4L) {
        stop(
          "The continuous phenotype does not contain enough ",
          "unique values to construct tertiles."
        )
      }
      
      plot_data[[phenotype_display_col]] <- cut(
        phenotype_values,
        breaks = breaks,
        include.lowest = TRUE,
        labels = c(
          "Lower tertile",
          "Middle tertile",
          "Upper tertile"
        )
      )
      
      display_type <- "categorical"
      
    } else {
      plot_data[[phenotype_display_col]] <-
        plot_data[[phenotype_col]]
      
      display_type <- "continuous"
    }
    
  } else {
    
    original_phenotype <- plot_data[[phenotype_col]]
    
    if (is.factor(original_phenotype)) {
      phenotype_levels <- levels(
        droplevels(original_phenotype)
      )
    } else {
      phenotype_levels <- unique(
        as.character(original_phenotype)
      )
    }
    
    plot_data[[phenotype_display_col]] <- factor(
      as.character(original_phenotype),
      levels = phenotype_levels
    )
    
    display_type <- "categorical"
  }
  
  # ============================================================
  # Optionally calculate group-specific CpG means
  # ============================================================
  
  if (plot_means) {
    
    plot_data <- plot_data |>
      dplyr::group_by(
        .data[[chr_col]],
        .data[[pos_col]],
        .data[[phenotype_display_col]]
      ) |>
      dplyr::summarise(
        mean_beta = mean(
          .data[[beta_col]],
          na.rm = TRUE
        ),
        mean_coverage = mean(
          .data[[coverage_col]],
          na.rm = TRUE
        ),
        median_coverage = stats::median(
          .data[[coverage_col]],
          na.rm = TRUE
        ),
        n_samples = sum(
          is.finite(.data[[beta_col]])
        ),
        beta_sd = dplyr::if_else(
          dplyr::n() > 1L,
          stats::sd(
            .data[[beta_col]],
            na.rm = TRUE
          ),
          NA_real_
        ),
        beta_se = dplyr::if_else(
          .data$n_samples > 1L,
          .data$beta_sd / sqrt(.data$n_samples),
          NA_real_
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        "{beta_col}" := .data$mean_beta,
        "{coverage_col}" := .data$mean_coverage
      )
    
    # Grouping strips some factor metadata; restore the original order.
    plot_data[[phenotype_display_col]] <- factor(
      as.character(plot_data[[phenotype_display_col]]),
      levels = phenotype_levels
    )
    
    # Group means already provide one profile per phenotype group.
    curve_method_for_plot <- "none"
    
  } else {
    curve_method_for_plot <- curve_method
  }
  
  # ============================================================
  # Prepare point-size variable
  # ============================================================
  
  if (
    !is.numeric(point_size_range) ||
    length(point_size_range) != 2L ||
    anyNA(point_size_range) ||
    point_size_range[[1L]] <= 0 ||
    point_size_range[[1L]] >= point_size_range[[2L]]
  ) {
    stop(
      "`point_size_range` must be a numeric vector of the form ",
      "c(minimum_size, maximum_size), for example c(1.5, 6)."
    )
  }
  
  if (
    !is.null(coverage_limits) &&
    (
      !is.numeric(coverage_limits) ||
      length(coverage_limits) != 2L ||
      anyNA(coverage_limits) ||
      coverage_limits[[1L]] >= coverage_limits[[2L]]
    )
  ) {
    stop(
      "`coverage_limits` must be NULL or a numeric vector ",
      "of the form c(lower, upper)."
    )
  }
  
  if (plot_means) {
    
    # Each point is a CpG x phenotype-group mean.
    # Its size represents mean sequencing coverage among the
    # samples contributing to that mean beta value.
    plot_data <- plot_data |>
      dplyr::mutate(
        .point_size_value = as.numeric(.data$mean_coverage)
      )
    
    size_legend_title <- "Mean coverage"
    
  } else {
    
    # Each point represents one observed sample x CpG measurement.
    # Its size represents that observation's sequencing coverage.
    plot_data <- plot_data |>
      dplyr::mutate(
        .point_size_value =
          as.numeric(.data[[coverage_col]])
      )
    
    size_legend_title <- "Coverage"
  }
  
  # Remove points without usable coverage from the size scale.
  plot_data <- plot_data |>
    dplyr::mutate(
      .point_size_value = dplyr::if_else(
        is.finite(.data$.point_size_value) &
          .data$.point_size_value >= 0,
        .data$.point_size_value,
        NA_real_
      )
    )
  
  # Optionally truncate extreme coverage values before mapping
  # them to point size.
  if (!is.null(coverage_limits)) {
    plot_data <- plot_data |>
      dplyr::mutate(
        .point_size_value = pmin(
          pmax(
            .data$.point_size_value,
            coverage_limits[[1L]]
          ),
          coverage_limits[[2L]]
        )
      )
    
    size_limits <- coverage_limits
    
  } else {
    size_limits <- NULL
  }
  
  # ============================================================
  # Create default categorical colors
  # ============================================================
  
  if (display_type == "categorical") {
    
    displayed_groups <- levels(
      droplevels(plot_data[[phenotype_display_col]])
    )
    
    if (is.null(group_colors)) {
      
      default_colors <- scales::hue_pal()(
        length(displayed_groups)
      )
      
      group_colors <- stats::setNames(
        default_colors,
        displayed_groups
      )
      
    } else {
      
      if (is.null(names(group_colors))) {
        if (length(group_colors) < length(displayed_groups)) {
          stop(
            "`group_colors` does not contain enough colors for ",
            length(displayed_groups),
            " phenotype groups."
          )
        }
        
        group_colors <- stats::setNames(
          group_colors[seq_along(displayed_groups)],
          displayed_groups
        )
      }
      
      missing_group_colors <- setdiff(
        displayed_groups,
        names(group_colors)
      )
      
      if (length(missing_group_colors) > 0L) {
        stop(
          "`group_colors` is missing colors for: ",
          paste(missing_group_colors, collapse = ", ")
        )
      }
      
      group_colors <- group_colors[displayed_groups]
    }
  }
  
  # ============================================================
  # Common plotting elements
  # ============================================================
  
  dmr_rectangle <- data.frame(
    xmin = dmr_start,
    xmax = dmr_end,
    ymin = -Inf,
    ymax = Inf
  )
  
  common_x_scale <- ggplot2::scale_x_continuous(
    limits = c(region_start, region_end),
    labels = scales::label_comma(),
    expand = ggplot2::expansion(mult = c(0.005, 0.005))
  )
  
  common_theme <- ggplot2::theme_classic(
    base_size = base_size
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        color = "black"
      )
    )
  
  # ============================================================
  # Methylation panel
  # ============================================================
  
  methylation_plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = dmr_rectangle,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      inherit.aes = FALSE,
      fill = dmr_highlight_color,
      alpha = dmr_highlight_alpha
    )
  
  if (display_type == "categorical") {
    
    if (plot_means && connect_means) {
      methylation_plot <- methylation_plot +
        ggplot2::geom_line(
          data = plot_data,
          ggplot2::aes(
            x = .data[[pos_col]],
            y = .data[[beta_col]],
            color = .data[[phenotype_display_col]],
            group = .data[[phenotype_display_col]]
          ),
          linewidth = mean_linewidth,
          alpha = 0.85,
          na.rm = TRUE
        )
    }
    
    methylation_plot <- methylation_plot +
      ggplot2::geom_point(
        data = plot_data,
        ggplot2::aes(
          x = .data[[pos_col]],
          y = .data[[beta_col]],
          color = .data[[phenotype_display_col]],
          size = .data$.point_size_value
        ),
        alpha = if (plot_means) 0.95 else point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        name = phenotype_col,
        values = group_colors,
        drop = FALSE
      )
    
    if (!plot_means && curve_method_for_plot == "loess") {
      methylation_plot <- methylation_plot +
        ggplot2::geom_smooth(
          data = plot_data,
          ggplot2::aes(
            x = .data[[pos_col]],
            y = .data[[beta_col]],
            color = .data[[phenotype_display_col]],
            group = .data[[phenotype_display_col]]
          ),
          method = "loess",
          formula = y ~ x,
          span = curve_span,
          se = FALSE,
          linewidth = curve_linewidth,
          alpha = curve_alpha,
          na.rm = TRUE
        )
    }
    
    if (!plot_means && curve_method_for_plot == "gam") {
      
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop(
          "Install `mgcv` to use `curve_method = 'gam'`."
        )
      }
      
      methylation_plot <- methylation_plot +
        ggplot2::geom_smooth(
          data = plot_data,
          ggplot2::aes(
            x = .data[[pos_col]],
            y = .data[[beta_col]],
            color = .data[[phenotype_display_col]],
            group = .data[[phenotype_display_col]]
          ),
          method = "gam",
          formula = y ~ s(x, bs = "cs"),
          se = FALSE,
          linewidth = curve_linewidth,
          alpha = curve_alpha,
          na.rm = TRUE
        )
    }
    
  } else {
    
    methylation_plot <- methylation_plot +
      ggplot2::geom_point(
        data = plot_data,
        ggplot2::aes(
          x = .data[[pos_col]],
          y = .data[[beta_col]],
          color = .data[[phenotype_display_col]],
          size = .data$.point_size_value
        ),
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_gradientn(
        name = phenotype_col,
        colors = continuous_palette
      )
    
    if (curve_method_for_plot == "loess") {
      methylation_plot <- methylation_plot +
        ggplot2::geom_smooth(
          data = plot_data,
          ggplot2::aes(
            x = .data[[pos_col]],
            y = .data[[beta_col]]
          ),
          method = "loess",
          formula = y ~ x,
          span = curve_span,
          se = FALSE,
          color = "black",
          linewidth = curve_linewidth,
          alpha = curve_alpha,
          na.rm = TRUE
        )
    }
    
    if (curve_method_for_plot == "gam") {
      
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop(
          "Install `mgcv` to use `curve_method = 'gam'`."
        )
      }
      
      methylation_plot <- methylation_plot +
        ggplot2::geom_smooth(
          data = plot_data,
          ggplot2::aes(
            x = .data[[pos_col]],
            y = .data[[beta_col]]
          ),
          method = "gam",
          formula = y ~ s(x, bs = "cs"),
          se = FALSE,
          color = "black",
          linewidth = curve_linewidth,
          alpha = curve_alpha,
          na.rm = TRUE
        )
    }
  }
  
  methylation_plot <- methylation_plot +
    ggplot2::scale_size_continuous(
      name = size_legend_title,
      range = point_size_range,
      limits = size_limits
    ) +
    common_x_scale +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(
        mult = c(0.01, 0.04)
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = if (plot_means) {
        "Mean DNA methylation"
      } else {
        "DNA methylation"
      },
      title = if (is.null(title)) {
        paste0(
          target_dmr_id,
          if (plot_means) " \u2014 group means" else ""
        )
      } else {
        title
      }
    ) +
    common_theme +
    ggplot2::theme(
      # Very light reference gridlines for CpG position and beta value.
      panel.grid.major = ggplot2::element_line(
        color = "grey92",
        linewidth = 0.30
      ),
      panel.grid.minor = ggplot2::element_line(
        color = "grey96",
        linewidth = 0.20
      ),
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.position = if (show_sample_legend) {
        "right"
      } else {
        "none"
      }
    )
  
  # ============================================================
  # Regulatory annotation track
  # ============================================================
  
  regulatory_plot <- NULL
  
  # Objects returned with the plot
  regulatory_region <- NULL
  cpg_region <- NULL
  genes_region <- NULL
  
  if (!is.null(regulatory_data)) {
    
    required_columns <- c(
      regulatory_chr_col,
      regulatory_start_col,
      regulatory_end_col,
      regulatory_type_col
    )
    
    missing_columns <- setdiff(
      required_columns,
      names(regulatory_data)
    )
    
    if (length(missing_columns) > 0L) {
      stop(
        "`regulatory_df` RDS is missing columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }
    
    regulatory_region <- regulatory_data |>
      dplyr::transmute(
        chr = as.character(
          .data[[regulatory_chr_col]]
        ),
        start = as.integer(
          .data[[regulatory_start_col]]
        ),
        end = as.integer(
          .data[[regulatory_end_col]]
        ),
        regulatory_type = tolower(
          as.character(
            .data[[regulatory_type_col]]
          )
        )
      ) |>
      dplyr::filter(
        .data$chr == .env$dmr_chr,
        .data$end >= .env$region_start,
        .data$start <= .env$region_end
      ) |>
      dplyr::mutate(
        start = pmax(.data$start, region_start),
        end = pmin(.data$end, region_end),
        y = 1
      )
    
    if (nrow(regulatory_region) > 0L) {
      
      regulatory_types <- unique(
        regulatory_region$regulatory_type
      )
      
      missing_colors <- setdiff(
        regulatory_types,
        names(regulatory_colors)
      )
      
      if (length(missing_colors) > 0L) {
        extra_colors <- scales::hue_pal()(
          length(missing_colors)
        )
        
        names(extra_colors) <- missing_colors
        
        regulatory_colors <- c(
          regulatory_colors,
          extra_colors
        )
      }
      
      regulatory_plot <- ggplot2::ggplot(
        regulatory_region
      ) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = .data$start,
            xend = .data$end,
            y = .data$y,
            yend = .data$y,
            color = .data$regulatory_type
          ),
          linewidth = 7,
          lineend = "butt"
        ) +
        ggplot2::scale_color_manual(
          name = "Regulatory element",
          values = regulatory_colors,
          breaks = regulatory_types,
          limits = regulatory_types,
          drop = TRUE
        ) +
        common_x_scale +
        ggplot2::scale_y_continuous(
          limits = c(0.5, 1.5)
        ) +
        ggplot2::labs(
          x = NULL,
          y = "Regulatory"
        ) +
        common_theme +
        ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.y = ggplot2::element_blank(),
          legend.position = "right"
        )
    }
  }
  
  # ============================================================
  # CpG-island track
  # ============================================================
  
  cpg_plot <- NULL
  
  if (!is.null(cpg_islands_data)) {
    
    required_columns <- c(
      cpg_chr_col,
      cpg_start_col,
      cpg_end_col
    )
    
    missing_columns <- setdiff(
      required_columns,
      names(cpg_islands_data)
    )
    
    if (length(missing_columns) > 0L) {
      stop(
        "`cpg_islands_df` RDS is missing columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }
    
    cpg_region <- cpg_islands_data |>
      dplyr::transmute(
        chr = as.character(
          .data[[cpg_chr_col]]
        ),
        start = as.integer(
          .data[[cpg_start_col]]
        ),
        end = as.integer(
          .data[[cpg_end_col]]
        ),
        label = if (
          !is.null(cpg_name_col) &&
          cpg_name_col %in% names(cpg_islands_data)
        ) {
          as.character(.data[[cpg_name_col]])
        } else {
          "CpG island"
        }
      ) |>
      dplyr::filter(
        .data$chr == .env$dmr_chr,
        .data$end >= .env$region_start,
        .data$start <= .env$region_end
      ) |>
      dplyr::mutate(
        start = pmax(.data$start, region_start),
        end = pmin(.data$end, region_end),
        y = 1
      )
    
    if (nrow(cpg_region) > 0L) {
      cpg_plot <- ggplot2::ggplot(cpg_region) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = .data$start,
            xend = .data$end,
            y = .data$y,
            yend = .data$y,
            color = "CpG island"
          ),
          linewidth = 7,
          lineend = "butt"
        ) +
        ggplot2::scale_color_manual(
          name = "CpG annotation",
          values = c(
            "CpG island" = cpg_island_color
          ),
          breaks = "CpG island",
          limits = "CpG island",
          drop = TRUE
        ) +
        common_x_scale +
        ggplot2::scale_y_continuous(
          limits = c(0.5, 1.5)
        ) +
        ggplot2::labs(
          x = NULL,
          y = "CpG island"
        ) +
        common_theme +
        ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.y = ggplot2::element_blank(),
          legend.position = "right"
        )
    }
  }
  
  # ============================================================
  # Gene-feature track
  # ============================================================
  
  gene_plot <- NULL
  
  if (!is.null(gene_features_data)) {
    
    required_columns <- c(
      gene_chr_col,
      gene_start_col,
      gene_end_col,
      gene_feature_col
    )
    
    missing_columns <- setdiff(
      required_columns,
      names(gene_features_data)
    )
    
    if (length(missing_columns) > 0L) {
      stop(
        "`gene_features_df` RDS is missing columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }
    
    if (!transcript_id_col %in% names(gene_features_data)) {
      stop(
        "`gene_features_df` must contain transcript column `",
        transcript_id_col, "`."
      )
    }
    
    has_gene_name <- gene_name_col %in% names(gene_features_data)
    has_gene_id <- gene_id_col %in% names(gene_features_data)
    has_strand <- gene_strand_col %in% names(gene_features_data)
    
    genes_region <- gene_features_data |>
      dplyr::mutate(
        .original_row = dplyr::row_number(),
        .plot_chr = as.character(
          .data[[gene_chr_col]]
        ),
        .plot_start = as.integer(
          .data[[gene_start_col]]
        ),
        .plot_end = as.integer(
          .data[[gene_end_col]]
        ),
        .plot_feature = tolower(
          as.character(.data[[gene_feature_col]])
        ),
        .plot_gene_name = if (has_gene_name) {
          as.character(.data[[gene_name_col]])
        } else if (has_gene_id) {
          as.character(.data[[gene_id_col]])
        } else {
          "Gene"
        },
        .plot_gene_id = if (has_gene_id) {
          as.character(.data[[gene_id_col]])
        } else {
          .data$.plot_gene_name
        },
        .plot_transcript_id = as.character(
          .data[[transcript_id_col]]
        ),
        .plot_strand = if (has_strand) {
          as.character(.data[[gene_strand_col]])
        } else {
          "*"
        }
      ) |>
      dplyr::filter(
        .data$.plot_chr == .env$dmr_chr,
        .data$.plot_end >= .env$region_start,
        .data$.plot_start <= .env$region_end
      )
    
    # ----------------------------------------------------------
    # Select one canonical transcript per gene
    # ----------------------------------------------------------
    
    if (nrow(genes_region) > 0L) {
      
      transcript_counts <- genes_region |>
        dplyr::distinct(
          .data$.plot_gene_id,
          .data$.plot_transcript_id
        ) |>
        dplyr::count(
          .data$.plot_gene_id,
          name = "n_transcripts"
        )
      
      genes_region <- genes_region |>
        dplyr::left_join(
          transcript_counts,
          by = ".plot_gene_id"
        )
      
      genes_with_multiple_transcripts <- genes_region |>
        dplyr::filter(.data$n_transcripts > 1L) |>
        dplyr::distinct(.data$.plot_gene_id) |>
        dplyr::pull(.data$.plot_gene_id)
      
      if (
        length(genes_with_multiple_transcripts) > 0L &&
        !canonical_transcript_id_col %in% names(genes_region)
      ) {
        stop(
          "Multiple transcripts occur for one or more genes in the ",
          "plotting window, but canonical transcript column `",
          canonical_transcript_id_col,
          "` is absent from `gene_features_df`."
        )
      }
      
      if (
        canonical_transcript_id_col %in%
        names(genes_region)
      ) {
        
        # Remove Ensembl version suffixes so, for example,
        # ENST00000274368 and ENST00000274368.9 match.
        genes_region <- genes_region |>
          dplyr::mutate(
            .plot_transcript_id_base = sub(
              "\\.[0-9]+$",
              "",
              .data$.plot_transcript_id
            ),
            .canonical_id_row = as.character(
              .data[[canonical_transcript_id_col]]
            ),
            .canonical_id_row = dplyr::na_if(
              trimws(.data$.canonical_id_row),
              ""
            ),
            .canonical_id_base_row = sub(
              "\\.[0-9]+$",
              "",
              .data$.canonical_id_row
            )
          )
        
        # canonical_transcript_id may be present only on the canonical
        # transcript's rows. Resolve it once per gene, then propagate it
        # to every row belonging to that gene.
        canonical_by_gene <- genes_region |>
          dplyr::filter(
            !is.na(.data$.canonical_id_base_row)
          ) |>
          dplyr::distinct(
            .data$.plot_gene_id,
            .data$.canonical_id_base_row
          ) |>
          dplyr::group_by(.data$.plot_gene_id) |>
          dplyr::summarise(
            .n_canonical_ids =
              dplyr::n_distinct(.data$.canonical_id_base_row),
            .resolved_canonical_id_base =
              dplyr::first(.data$.canonical_id_base_row),
            .groups = "drop"
          )
        
        ambiguous_canonical_genes <- canonical_by_gene |>
          dplyr::filter(.data$.n_canonical_ids > 1L) |>
          dplyr::pull(.data$.plot_gene_id)
        
        if (length(ambiguous_canonical_genes) > 0L) {
          stop(
            "More than one canonical transcript ID was found for: ",
            paste(ambiguous_canonical_genes, collapse = ", ")
          )
        }
        
        genes_region <- genes_region |>
          dplyr::left_join(
            canonical_by_gene |>
              dplyr::select(
                .data$.plot_gene_id,
                .data$.resolved_canonical_id_base
              ),
            by = ".plot_gene_id"
          )
        
        # Check once per gene, rather than checking every transcript row.
        missing_canonical_genes <- genes_region |>
          dplyr::filter(.data$n_transcripts > 1L) |>
          dplyr::distinct(
            .data$.plot_gene_id,
            .data$.resolved_canonical_id_base
          ) |>
          dplyr::filter(
            is.na(.data$.resolved_canonical_id_base) |
              !nzchar(.data$.resolved_canonical_id_base)
          ) |>
          dplyr::pull(.data$.plot_gene_id)
        
        if (length(missing_canonical_genes) > 0L) {
          stop(
            "The following genes have multiple transcripts but no ",
            "usable canonical transcript ID: ",
            paste(missing_canonical_genes, collapse = ", ")
          )
        }
        
        unmatched_canonical_genes <- genes_region |>
          dplyr::filter(.data$n_transcripts > 1L) |>
          dplyr::group_by(.data$.plot_gene_id) |>
          dplyr::summarise(
            canonical_found = any(
              .data$.plot_transcript_id_base ==
                .data$.resolved_canonical_id_base,
              na.rm = TRUE
            ),
            .groups = "drop"
          ) |>
          dplyr::filter(!.data$canonical_found) |>
          dplyr::pull(.data$.plot_gene_id)
        
        if (length(unmatched_canonical_genes) > 0L) {
          stop(
            "The canonical transcript ID did not match any transcript ",
            "rows for: ",
            paste(unmatched_canonical_genes, collapse = ", ")
          )
        }
        
        # For genes with multiple transcripts, retain only the canonical
        # transcript. A single-transcript gene is retained as-is.
        genes_region <- genes_region |>
          dplyr::mutate(
            .is_canonical = dplyr::if_else(
              .data$n_transcripts == 1L,
              TRUE,
              .data$.plot_transcript_id_base ==
                .data$.resolved_canonical_id_base,
              missing = FALSE
            )
          ) |>
          dplyr::filter(.data$.is_canonical)
      }
    }
    
    # Clip annotations to the plotting window only after canonical
    # transcript selection.
    genes_region <- genes_region |>
      dplyr::mutate(
        .plot_start_original = .data$.plot_start,
        .plot_end_original = .data$.plot_end,
        .plot_start = pmax(
          .data$.plot_start,
          .env$region_start
        ),
        .plot_end = pmin(
          .data$.plot_end,
          .env$region_end
        ),
        .plot_strand = dplyr::case_when(
          .data$.plot_strand %in% c("+", "1", 1) ~ "+",
          .data$.plot_strand %in% c("-", "-1", -1) ~ "-",
          TRUE ~ "*"
        ),
        .feature_class = dplyr::case_when(
          grepl(
            "tss|transcription.start",
            .data$.plot_feature,
            ignore.case = TRUE
          ) ~ "TSS",
          grepl(
            "exon",
            .data$.plot_feature,
            ignore.case = TRUE
          ) ~ "Exon",
          grepl(
            "intron",
            .data$.plot_feature,
            ignore.case = TRUE
          ) ~ "Intron",
          TRUE ~ "Other"
        ),
        .transcript_label = paste0(
          .data$.plot_gene_name,
          " | ",
          .data$.plot_transcript_id
        )
      )
    
    if (nrow(genes_region) > 0L) {
      
      transcript_order <- genes_region |>
        dplyr::group_by(
          .data$.transcript_label,
          .data$.plot_gene_name
        ) |>
        dplyr::summarise(
          .transcript_start = min(
            .data$.plot_start,
            na.rm = TRUE
          ),
          .groups = "drop"
        ) |>
        dplyr::arrange(
          .data$.transcript_start,
          .data$.plot_gene_name,
          .data$.transcript_label
        ) |>
        dplyr::slice_head(
          n = max_gene_transcripts
        ) |>
        dplyr::pull(.data$.transcript_label) |>
        unique()
      
      genes_region <- genes_region |>
        dplyr::filter(
          .data$.transcript_label %in%
            .env$transcript_order
        ) |>
        dplyr::mutate(
          .transcript_label = factor(
            .data$.transcript_label,
            levels = rev(unique(.env$transcript_order))
          ),
          .transcript_y = as.numeric(.data$.transcript_label)
        )
      
      # ----------------------------------------------------------
      # Generate TSS arrows
      # ----------------------------------------------------------
      
      plotted_window_width <- region_end - region_start + 1L
      
      tss_arrow_length <- max(
        1L,
        round(
          plotted_window_width *
            tss_arrow_length_prop
        )
      )
      
      tss_arrows <- genes_region |>
        dplyr::filter(.data$.feature_class == "TSS") |>
        dplyr::mutate(
          # Use the biologically relevant TSS coordinate before clipping.
          .tss_position = dplyr::if_else(
            .data$.plot_strand == "-",
            .data$.plot_end_original,
            .data$.plot_start_original
          ),
          .arrow_start = .data$.tss_position,
          .arrow_end = dplyr::case_when(
            .data$.plot_strand == "+" ~
              .data$.tss_position + tss_arrow_length,
            .data$.plot_strand == "-" ~
              .data$.tss_position - tss_arrow_length,
            TRUE ~ .data$.tss_position
          ),
          .arrow_start = pmin(
            pmax(.data$.arrow_start, region_start),
            region_end
          ),
          .arrow_end = pmin(
            pmax(.data$.arrow_end, region_start),
            region_end
          )
        ) |>
        dplyr::filter(
          .data$.tss_position >= .env$region_start,
          .data$.tss_position <= .env$region_end
        ) |>
        dplyr::distinct(
          .data$.transcript_label,
          .data$.plot_strand,
          .data$.tss_position,
          .data$.arrow_start,
          .data$.arrow_end
        )
      
      tss_arrows <- tss_arrows |>
        dplyr::mutate(
          .transcript_y = as.numeric(.data$.transcript_label),
          .tss_arrow_y = .data$.transcript_y +
            tss_vertical_nudge
        )
      # ----------------------------------------------------------
      # Generate equally spaced intron arrows
      # ----------------------------------------------------------
      
      if (is.null(intron_arrow_spacing)) {
        resolved_intron_arrow_spacing <- max(
          25L,
          round(plotted_window_width / 30)
        )
      } else {
        if (
          length(intron_arrow_spacing) != 1L ||
          !is.numeric(intron_arrow_spacing) ||
          is.na(intron_arrow_spacing) ||
          intron_arrow_spacing <= 0
        ) {
          stop(
            "`intron_arrow_spacing` must be NULL or one positive ",
            "number representing base pairs."
          )
        }
        
        resolved_intron_arrow_spacing <- as.numeric(
          intron_arrow_spacing
        )
      }
      
      intron_arrow_length <- max(
        1,
        resolved_intron_arrow_spacing *
          intron_arrow_length_prop
      )
      
      make_intron_arrows <- function(
    start,
    end,
    strand,
    transcript_label
      ) {
        
        segment_width <- end - start
        
        if (
          !is.finite(segment_width) ||
          segment_width <
          resolved_intron_arrow_spacing ||
          !strand %in% c("+", "-")
        ) {
          return(NULL)
        }
        
        arrow_centers <- seq(
          from = start +
            resolved_intron_arrow_spacing / 2,
          to = end -
            resolved_intron_arrow_spacing / 2,
          by = resolved_intron_arrow_spacing
        )
        
        if (length(arrow_centers) == 0L) {
          return(NULL)
        }
        
        if (strand == "+") {
          arrow_start <- arrow_centers -
            intron_arrow_length / 2
          arrow_end <- arrow_centers +
            intron_arrow_length / 2
        } else {
          arrow_start <- arrow_centers +
            intron_arrow_length / 2
          arrow_end <- arrow_centers -
            intron_arrow_length / 2
        }
        
        data.frame(
          .transcript_label = transcript_label,
          .plot_strand = strand,
          .arrow_start = arrow_start,
          .arrow_end = arrow_end,
          stringsAsFactors = FALSE
        )
      }
      
      intron_rows <- genes_region |>
        dplyr::filter(
          .data$.feature_class %in%
            c("Intron", "Other")
        )
      
      if (nrow(intron_rows) > 0L) {
        intron_arrows <- lapply(
          seq_len(nrow(intron_rows)),
          function(i) {
            make_intron_arrows(
              start = intron_rows$.plot_start[[i]],
              end = intron_rows$.plot_end[[i]],
              strand = intron_rows$.plot_strand[[i]],
              transcript_label = as.character(
                intron_rows$.transcript_label[[i]]
              )
            )
          }
        )
        
        intron_arrows <- dplyr::bind_rows(intron_arrows)
        
        if (nrow(intron_arrows) > 0L) {
          intron_arrows$.transcript_label <- factor(
            intron_arrows$.transcript_label,
            levels = levels(
              genes_region$.transcript_label
            )
          )
        }
      } else {
        intron_arrows <- data.frame()
      }
      intron_arrows <- intron_arrows |>
        dplyr::mutate(
          .transcript_y = as.numeric(.data$.transcript_label)
        )
      
      gene_features_present <- c(
        if (any(
          genes_region$.feature_class == "Exon",
          na.rm = TRUE
        )) {
          "Exon"
        },
        if (any(
          genes_region$.feature_class %in% c("Intron", "Other"),
          na.rm = TRUE
        )) {
          "Intron"
        },
        if (nrow(tss_arrows) > 0L) {
          "TSS"
        }
      ) |>
        unique()
      
      gene_feature_colors <- c(
        "Exon" = exon_color,
        "Intron" = intron_color,
        "TSS" = tss_color
      )
      
      gene_plot <- ggplot2::ggplot() +
        
        # Baseline gene/intron segments
        ggplot2::geom_segment(
          data = genes_region |>
            dplyr::filter(
              .data$.feature_class %in% c("Intron", "Other")
            ) |>
            dplyr::mutate(
              .legend_feature = "Intron"
            ),
          ggplot2::aes(
            x = .data$.plot_start,
            xend = .data$.plot_end,
            y = .data$.transcript_y,
            yend = .data$.transcript_y,
            color = .data$.legend_feature
          ),
          linewidth = 0.65,
          lineend = "butt"
        ) +
        
        # Exon blocks
        ggplot2::geom_segment(
          data = genes_region |>
            dplyr::filter(
              .data$.feature_class == "Exon"
            ),
          ggplot2::aes(
            x = .data$.plot_start,
            xend = .data$.plot_end,
            y = .data$.transcript_y,
            yend = .data$.transcript_y,
            color = .data$.feature_class
          ),
          linewidth = 4,
          lineend = "butt"
        )
      
      # Add directional arrows along plotted introns.
      if (nrow(intron_arrows) > 0L) {
        gene_plot <- gene_plot +
          ggplot2::geom_segment(
            data = intron_arrows |>
              dplyr::mutate(
                .legend_feature = "Intron"
              ),
            ggplot2::aes(
              x = .data$.arrow_start,
              xend = .data$.arrow_end,
              y = .data$.transcript_y,
              yend = .data$.transcript_y,
              color = .data$.legend_feature
            ),
            linewidth = 0.55,
            arrow = grid::arrow(
              type = "closed",
              ends = "last",
              length = grid::unit(
                annotation_arrow_size,
                "cm"
              )
            )
          )
      }
      
      # Add a TSS line and directional arrow for each transcript.
      # Add an L-shaped TSS marker:
      # vertical stem followed by a horizontal directional arrow.
      if (nrow(tss_arrows) > 0L) {
        
        directional_tss_arrows <- tss_arrows |>
          dplyr::filter(
            .data$.plot_strand %in% c("+", "-"),
            .data$.arrow_start != .data$.arrow_end
          )
        
        gene_plot <- gene_plot +
          
          # Vertical stem rising from the transcript lane.
          ggplot2::geom_segment(
            data = tss_arrows |>
              dplyr::mutate(
                .legend_feature = "TSS"
              ),
            ggplot2::aes(
              x = .data$.tss_position,
              xend = .data$.tss_position,
              y = .data$.transcript_y,
              yend = .data$.tss_arrow_y,
              color = .data$.legend_feature
            ),
            linewidth = 0.85,
            lineend = "butt"
          )
        
        if (nrow(directional_tss_arrows) > 0L) {
          
          gene_plot <- gene_plot +
            
            # Horizontal arm beginning at the top of the vertical stem.
            # The arrowhead occurs only at the terminal end.
            ggplot2::geom_segment(
              data = directional_tss_arrows |>
                dplyr::mutate(
                  .legend_feature = "TSS"
                ),
              ggplot2::aes(
                x = .data$.arrow_start,
                xend = .data$.arrow_end,
                y = .data$.tss_arrow_y,
                yend = .data$.tss_arrow_y,
                color = .data$.legend_feature
              ),
              linewidth = 0.85,
              lineend = "butt",
              arrow = grid::arrow(
                type = "closed",
                ends = "last",
                length = grid::unit(
                  annotation_arrow_size,
                  "cm"
                )
              )
            )
        }
      }
      
      gene_plot <- gene_plot +
        ggplot2::scale_color_manual(
          name = "Gene annotation",
          values = gene_feature_colors,
          breaks = gene_features_present,
          limits = gene_features_present,
          drop = TRUE,
          guide = ggplot2::guide_legend(
            override.aes = list(
              linewidth = unname(
                c(
                  "Exon" = 4,
                  "Intron" = 1,
                  "TSS" = 1
                )[gene_features_present]
              )
            )
          )
        ) +
        common_x_scale +
        ggplot2::labs(
          x = NULL,
          y = NULL
        ) +
        ggplot2::scale_y_continuous(
          breaks = seq_along(
            levels(genes_region$.transcript_label)
          ),
          labels = levels(
            genes_region$.transcript_label
          ),
          expand = ggplot2::expansion(
            add = c(0.25, tss_vertical_nudge + 0.25)
          )
        ) +
        common_theme +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(
            size = base_size * 0.70
          ),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.y = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          legend.position = "right"
        )
    }
  }
  
  # ============================================================
  # Combine available tracks
  # ============================================================
  
  plot_list <- list(
    methylation = methylation_plot
  )
  
  height_list <- c(
    methylation = relative_track_heights[["methylation"]]
  )
  
  if (!is.null(gene_plot)) {
    plot_list$genes <- gene_plot
    height_list[["genes"]] <-
      relative_track_heights[["genes"]]
  }
  
  if (!is.null(cpg_plot)) {
    plot_list$cpg <- cpg_plot
    height_list[["cpg"]] <-
      relative_track_heights[["cpg"]]
  }
  
  if (!is.null(regulatory_plot)) {
    plot_list$regulatory <- regulatory_plot
    height_list[["regulatory"]] <-
      relative_track_heights[["regulatory"]]
  }
  
  combined_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = 1,
    heights = unname(height_list),
    guides = "collect"
  ) &
    ggplot2::theme(
      legend.position = if (show_sample_legend) {
        "right"
      } else {
        "none"
      }
    )
  
  # Add the genomic-position label to the bottom panel.
  bottom_panel_name <- names(plot_list)[length(plot_list)]
  
  plot_list[[bottom_panel_name]] <-
    plot_list[[bottom_panel_name]] +
    ggplot2::labs(
      x = paste0(
        dmr_chr,
        " position (bp)"
      )
    )
  
  combined_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = 1,
    heights = unname(height_list),
    guides = "collect"
  ) &
    ggplot2::theme(
      legend.position = if (show_sample_legend) {
        "right"
      } else {
        "none"
      }
    )
  
  # ============================================================
  # Return plot and all plotted data
  # ============================================================
  
  plot_region <- data.frame(
    chr = dmr_chr,
    start = region_start,
    end = region_end,
    stringsAsFactors = FALSE
  )
  
  # Return zero-row data frames instead of NULL when an annotation
  # file was supplied but no feature overlapped the window.
  if (!is.null(regulatory_data) && is.null(regulatory_region)) {
    regulatory_region <- regulatory_data[0, , drop = FALSE]
  }
  
  if (!is.null(cpg_islands_data) && is.null(cpg_region)) {
    cpg_region <- cpg_islands_data[0, , drop = FALSE]
  }
  
  if (!is.null(gene_features_data) && is.null(genes_region)) {
    genes_region <- gene_features_data[0, , drop = FALSE]
  }
  
  return(
    list(
      plot = combined_plot,
      methylation_plot_data = plot_data,
      regulatory_annotations = regulatory_region,
      cpg_annotations = cpg_region,
      gene_annotations = genes_region,
      selected_dmr = selected_dmr,
      plot_region = plot_region
    )
  )
}
