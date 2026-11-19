.CPGDB_DEFAULTS <- list(
  db_file = NULL,
  table_name = "cpg_metadata_master",
  lib_path = NULL
)

cpgdb_defaults <- function() {
  list(
    db_file    = getOption("cpgdb.db_file", Sys.getenv("CPGDB_FILE", unset = "")),
    table_name = getOption("cpgdb.table_name", .CPGDB_DEFAULTS$table_name),
    lib_path   = getOption("cpgdb.lib_path", .CPGDB_DEFAULTS$lib_path)
  )
}

.REMOVE_RULES <- c(
  xy              = "remove_xy",
  all_dbsnp       = "remove_all_dbsnp",
  ct_ga           = "remove_ct_ga",
  maf1            = "remove_maf1",
  maf5            = "remove_maf5",
  all_dbsnp_no_xy = "remove_all_dbsnp_no_xy",
  ct_ga_no_xy     = "remove_ct_ga_no_xy",
  maf1_no_xy      = "remove_maf1_no_xy",
  maf5_no_xy      = "remove_maf5_no_xy"
)

.COLUMN_PREFIXES <- c("", "anno_", "var_", "arraymap_")
.MAX_SQL_PARAMS  <- 30000L   # SQLite allows 32766 bound parameters per statement

# ---- Internal helpers ---------------------------------------------------------
# DBI/RSQLite come from the user's own library; the lab library is only appended
# (never prepended) as a fallback, so sourcing this file doesn't change which
# versions of other packages a user's session loads.
.cpgdb_runtime <- function() {
  need <- c("DBI", "RSQLite")
  have <- function() vapply(need, requireNamespace, logical(1), quietly = TRUE)
  ok <- have()
  if (!all(ok)) {
    lib <- cpgdb_defaults()$lib_path
    if (!is.null(lib) && dir.exists(lib)) {
      .libPaths(unique(c(.libPaths(), lib)))
      ok <- have()
    }
  }
  if (!all(ok)) stop("Missing package(s): ", paste(need[!ok], collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

.cpgdb_path <- function(db_file = NULL) {
  db_file <- db_file %||% cpgdb_defaults()$db_file
  if (length(db_file) != 1L || is.na(db_file) || !nzchar(db_file) || !file.exists(db_file)) {
    stop("CpG metadata DB not found: ", db_file,
         "\nSet options(cpgdb.db_file = '...') or the CPGDB_FILE environment variable.", call. = FALSE)
  }
  db_file
}

#' Open a read-only connection. Close it with DBI::dbDisconnect(con).
cpgdb_connect <- function(db_file = NULL) {
  .cpgdb_runtime()
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = .cpgdb_path(db_file),
                        flags = RSQLite::SQLITE_RO, synchronous = NULL)
  DBI::dbExecute(con, "PRAGMA query_only = ON;")
  DBI::dbExecute(con, "PRAGMA cache_size = -262144;")    # 256 MB page cache
  DBI::dbExecute(con, "PRAGMA mmap_size = 1073741824;")  # memory-map up to 1 GB
  con
}

.with_con <- function(fn, db_file = NULL) {
  con <- cpgdb_connect(db_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  fn(con)
}

.check_table <- function(con, table_name) {
  if (!DBI::dbExistsTable(con, table_name)) {
    stop("Table '", table_name, "' not found. Available: ",
         paste(DBI::dbListTables(con), collapse = ", "), call. = FALSE)
  }
  DBI::dbListFields(con, table_name)
}

# Map short/full column names to the actual names in `fields`.
.resolve_columns <- function(columns, fields) {
  vapply(columns, function(col) {
    hit <- intersect(paste0(.COLUMN_PREFIXES, col), fields)
    if (length(hit) == 0) {
      stop("Column '", col, "' not found. See get_cpg_column_names().", call. = FALSE)
    }
    hit[1]
  }, character(1), USE.NAMES = FALSE)
}

# First existing column among candidates (or NA).
.first_col <- function(candidates, fields) intersect(candidates, fields)[1]

.ids_from_matrix <- function(mat) {
  if (!is.matrix(mat) && !is.data.frame(mat)) stop("'mat' must be a matrix or data.frame.", call. = FALSE)
  ids <- rownames(mat)
  if (is.null(ids) || (is.data.frame(mat) && .row_names_info(mat) < 0)) {
    stop("'mat' must have rownames that are CpG IDs (chr.pos.end).", call. = FALSE)
  }
  ids
}

# TRUE where any ';'-separated token of x is in `set` (vectorised).
.any_token_in <- function(x, set) {
  out <- logical(length(x))
  ok <- !is.na(x) & nzchar(x)
  if (!any(ok)) return(out)
  parts <- strsplit(x[ok], ";", fixed = TRUE)
  hits <- rep(which(ok), lengths(parts))[unlist(parts, use.names = FALSE) %in% set]
  out[unique(hits)] <- TRUE
  out
}

# ---- 1) Core fetch ----------------------------------------------------------------
#' Metadata for the requested CpGs, one row per input ID in input order
#' (duplicates kept). Adds query_order and requested_cpg_id columns.
#' @param columns  columns to return (full or short names); NULL = all.
#' @param strict   TRUE: error if any ID is missing; FALSE: warn and return NA rows.
fetch_cpg_meta <- function(cpg_ids,
                           columns = NULL,
                           strict = TRUE,
                           db_file = NULL,
                           table_name = NULL) {
  cpg_ids <- as.character(cpg_ids)
  if (anyNA(cpg_ids)) stop("'cpg_ids' contains NA values.", call. = FALSE)
  table_name <- table_name %||% cpgdb_defaults()$table_name
  
  .with_con(function(con) {
    fields <- .check_table(con, table_name)
    if (!"cpg_id" %in% fields) stop("Table '", table_name, "' has no cpg_id column.", call. = FALSE)
    cols <- if (is.null(columns)) fields else unique(c("cpg_id", .resolve_columns(columns, fields)))
    sel <- paste(DBI::dbQuoteIdentifier(con, cols), collapse = ", ")
    tq  <- DBI::dbQuoteIdentifier(con, table_name)
    
    uniq <- unique(cpg_ids)
    meta <- DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s LIMIT 0", sel, tq))  # typed empty frame
    if (length(uniq) > 0) {
      chunks <- split(uniq, ceiling(seq_along(uniq) / .MAX_SQL_PARAMS))
      meta <- do.call(rbind, lapply(chunks, function(ids) {
        DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s WHERE cpg_id IN (%s)", sel, tq,
                                     paste(rep("?", length(ids)), collapse = ",")),
                        params = as.list(ids))
      }))
    }
    
    m <- match(cpg_ids, meta$cpg_id)
    missing <- unique(cpg_ids[is.na(m)])
    if (length(missing)) {
      txt <- sprintf(
        "%d of %d unique CpG IDs not found in '%s' (e.g. %s). They are off-target or not in chr.pos.end format.",
        length(missing), length(uniq), table_name, paste(utils::head(missing, 3), collapse = ", "))
      if (strict) stop(txt, call. = FALSE) else warning(txt, call. = FALSE)
    }
    
    out <- cbind(data.frame(query_order = seq_along(cpg_ids), requested_cpg_id = cpg_ids,
                            stringsAsFactors = FALSE),
                 meta[m, , drop = FALSE])
    rownames(out) <- NULL
    out
  }, db_file = db_file)
}

# ---- 2) Annotate CpGs -------------------------------------------------------------
#' Like fetch_cpg_meta(), without query_order; rownames are the CpG IDs when unique.
annotate_cpgs <- function(cpg_ids,
                          columns = NULL,
                          strict = TRUE,
                          db_file = NULL,
                          table_name = NULL) {
  out <- fetch_cpg_meta(cpg_ids, columns = columns, strict = strict,
                        db_file = db_file, table_name = table_name)
  out$query_order <- NULL
  if (!anyDuplicated(out$requested_cpg_id)) rownames(out) <- out$requested_cpg_id
  out
}

#' Metadata for the rownames of a beta/M-value matrix, in row order.
get_metadata_for_beta <- function(mat,
                                  columns = NULL,
                                  strict = TRUE,
                                  db_file = NULL,
                                  table_name = NULL) {
  annotate_cpgs(.ids_from_matrix(mat), columns = columns, strict = strict,
                db_file = db_file, table_name = table_name)
}

# ---- 3) Rule-based variant filters ------------------------------------------------
#' Rule short name -> base column name (e.g. "ct_ga_no_xy" -> "remove_ct_ga_no_xy").
resolve_remove_rule <- function(rule = names(.REMOVE_RULES)) {
  .REMOVE_RULES[[match.arg(rule)]]
}

#' Named logical vector: TRUE = remove. `rule` is a short name (see .REMOVE_RULES)
#' or any boolean column name. With strict = FALSE, CpGs not in the DB are flagged
#' for removal (they have no variant QC).
get_remove_vector <- function(cpg_ids,
                              rule = "ct_ga_no_xy",
                              strict = TRUE,
                              db_file = NULL,
                              table_name = NULL) {
  col <- if (rule %in% names(.REMOVE_RULES)) .REMOVE_RULES[[rule]] else rule
  ann <- fetch_cpg_meta(cpg_ids, columns = col, strict = strict,
                        db_file = db_file, table_name = table_name)
  value_col <- setdiff(colnames(ann), c("query_order", "requested_cpg_id", "cpg_id"))
  if (length(value_col) == 0) value_col <- "cpg_id"   # rule column was cpg_id itself (unlikely)
  remove <- as.logical(ann[[value_col[1]]])
  remove[is.na(remove)] <- TRUE
  names(remove) <- ann$requested_cpg_id
  remove
}

#' Drop matrix rows flagged by a remove rule.
filter_cpg_matrix <- function(mat,
                              rule = "ct_ga_no_xy",
                              strict = TRUE,
                              verbose = TRUE,
                              db_file = NULL,
                              table_name = NULL) {
  remove <- get_remove_vector(.ids_from_matrix(mat), rule = rule, strict = strict,
                              db_file = db_file, table_name = table_name)
  out <- mat[!remove, , drop = FALSE]
  if (verbose) message(sprintf("Applied rule '%s': retained %d of %d rows.", rule, nrow(out), nrow(mat)))
  out
}

# ---- 4) Flexible metadata filter --------------------------------------------------
#' Keep matrix rows whose metadata pass every filter given (filters left NULL are ignored).
#' @param maf_pct_max keep CpGs with 1000G MAF <= this (in percent). CpGs with no
#'                    1000G SNV count as MAF 0 and are kept.
#' @param maf_pct_min keep CpGs with a 1000G SNV of MAF >= this (in percent).
#' @param has_snp     TRUE/FALSE: dbSNP155 overlap at the CpG.
#' @param cpg_context_include/_exclude  CpG_island, CpG_shore, CpG_shelf, OpenSea.
#' @param epigenetic_include/_exclude   epigenetic_primary values.
#' @param genic_include/_exclude        genic_primary values (promoter, 5UTR, exon, ...).
#' @param gene_symbols_include/_exclude matched against genic_gene_symbols.
#' @param tss_abs_max/_min              absolute distance (bp) to nearest TSS.
#' @param return  "matrix", "metadata" or "both".
subset_matrix_by_metadata <- function(mat,
                                      strict = TRUE,
                                      verbose = TRUE,
                                      maf_pct_max = NULL,
                                      maf_pct_min = NULL,
                                      has_snp = NULL,
                                      cpg_context_include = NULL,
                                      cpg_context_exclude = NULL,
                                      epigenetic_include = NULL,
                                      epigenetic_exclude = NULL,
                                      genic_include = NULL,
                                      genic_exclude = NULL,
                                      gene_symbols_include = NULL,
                                      gene_symbols_exclude = NULL,
                                      tss_abs_max = NULL,
                                      tss_abs_min = NULL,
                                      return = c("matrix", "metadata", "both"),
                                      db_file = NULL,
                                      table_name = NULL) {
  return <- match.arg(return)
  ids <- .ids_from_matrix(mat)
  
  # Columns each filter needs (short names, resolved against the table's prefixes).
  need <- c(
    if (!is.null(maf_pct_max) || !is.null(maf_pct_min)) "maf_pct",
    if (!is.null(has_snp)) "all_dbsnp",
    if (!is.null(cpg_context_include) || !is.null(cpg_context_exclude)) "cpg_context",
    if (!is.null(epigenetic_include) || !is.null(epigenetic_exclude)) "epigenetic_primary",
    if (!is.null(genic_include) || !is.null(genic_exclude)) "genic_primary",
    if (!is.null(gene_symbols_include) || !is.null(gene_symbols_exclude)) "genic_gene_symbols",
    if (!is.null(tss_abs_max) || !is.null(tss_abs_min)) "tss_distance_bp_abs"
  )
  # Matrix-only output needs just the filter columns; metadata output gets everything.
  ann <- fetch_cpg_meta(ids, columns = if (return == "matrix") need else NULL, strict = strict,
                        db_file = db_file, table_name = table_name)
  col <- function(short) ann[[.resolve_columns(short, colnames(ann))]]
  
  keep <- !is.na(ann$cpg_id)   # CpGs missing from the DB never pass (strict = FALSE)
  
  if (!is.null(maf_pct_max)) { maf <- col("maf_pct"); keep <- keep & (is.na(maf) | maf <= maf_pct_max) }
  if (!is.null(maf_pct_min)) { maf <- col("maf_pct"); keep <- keep & !is.na(maf) & maf >= maf_pct_min }
  if (!is.null(has_snp))     keep <- keep & (as.logical(col("all_dbsnp")) %in% has_snp)
  if (!is.null(cpg_context_include)) keep <- keep & col("cpg_context") %in% cpg_context_include
  if (!is.null(cpg_context_exclude)) keep <- keep & !(col("cpg_context") %in% cpg_context_exclude)
  if (!is.null(epigenetic_include))  keep <- keep & col("epigenetic_primary") %in% epigenetic_include
  if (!is.null(epigenetic_exclude))  keep <- keep & !(col("epigenetic_primary") %in% epigenetic_exclude)
  if (!is.null(genic_include))       keep <- keep & col("genic_primary") %in% genic_include
  if (!is.null(genic_exclude))       keep <- keep & !(col("genic_primary") %in% genic_exclude)
  if (!is.null(gene_symbols_include)) keep <- keep & .any_token_in(col("genic_gene_symbols"), gene_symbols_include)
  if (!is.null(gene_symbols_exclude)) keep <- keep & !.any_token_in(col("genic_gene_symbols"), gene_symbols_exclude)
  if (!is.null(tss_abs_max)) { d <- col("tss_distance_bp_abs"); keep <- keep & !is.na(d) & d <= tss_abs_max }
  if (!is.null(tss_abs_min)) { d <- col("tss_distance_bp_abs"); keep <- keep & !is.na(d) & d >= tss_abs_min }
  
  keep[is.na(keep)] <- FALSE
  mat_out <- mat[keep, , drop = FALSE]          # by position: safe with duplicated rownames
  ann_out <- ann[keep, setdiff(colnames(ann), "query_order"), drop = FALSE]
  if (!anyDuplicated(ann_out$requested_cpg_id)) rownames(ann_out) <- ann_out$requested_cpg_id
  
  if (verbose) message(sprintf("Metadata filter retained %d of %d CpGs.", sum(keep), length(keep)))
  switch(return,
         matrix   = mat_out,
         metadata = ann_out,
         both     = list(matrix = mat_out, metadata = ann_out))
}

# ---- 5) Region query ------------------------------------------------------------------
#' CpGs in one or more hg38 regions (1-based, inclusive), using the position index.
#' Give regions as strings ("chr7:27,100,000-27,250,000") or as chrom/start/end vectors.
#' With several regions a region_id column says which region each row came from.
get_cpgs_in_region <- function(chrom,
                               start = NULL,
                               end = NULL,
                               columns = NULL,
                               db_file = NULL,
                               table_name = NULL) {
  if (is.null(start) || is.null(end)) {
    reg <- gsub(",", "", as.character(chrom))
    ok <- grepl("^[^:]+:[0-9]+-[0-9]+$", reg)
    if (!all(ok)) stop("Regions must look like 'chr7:100-200'; bad: ", paste(reg[!ok], collapse = ", "), call. = FALSE)
    chrom <- sub(":.*$", "", reg)
    start <- as.integer(sub("^[^:]+:([0-9]+)-.*$", "\\1", reg))
    end   <- as.integer(sub("^.*-([0-9]+)$", "\\1", reg))
  }
  n <- max(length(chrom), length(start), length(end))
  chrom <- rep_len(as.character(chrom), n)
  start <- rep_len(as.integer(start), n)
  end   <- rep_len(as.integer(end), n)
  if (anyNA(start) || anyNA(end) || any(end < start)) stop("Invalid start/end.", call. = FALSE)
  table_name <- table_name %||% cpgdb_defaults()$table_name
  
  .with_con(function(con) {
    fields <- .check_table(con, table_name)
    chr_col <- .first_col(c("anno_chrom", "chrom", "chr", "var_chr"), fields)
    pos_col <- .first_col(c("anno_pos", "pos", "var_pos"), fields)
    if (is.na(chr_col) || is.na(pos_col)) stop("Table '", table_name, "' has no chrom/pos columns.", call. = FALSE)
    cols <- if (is.null(columns)) fields else unique(c("cpg_id", .resolve_columns(columns, fields)))
    sql <- sprintf("SELECT %s FROM %s WHERE %s = ? AND %s BETWEEN ? AND ? ORDER BY %s",
                   paste(DBI::dbQuoteIdentifier(con, cols), collapse = ", "),
                   DBI::dbQuoteIdentifier(con, table_name),
                   DBI::dbQuoteIdentifier(con, chr_col), DBI::dbQuoteIdentifier(con, pos_col),
                   DBI::dbQuoteIdentifier(con, pos_col))
    res <- lapply(seq_len(n), function(i) {
      r <- DBI::dbGetQuery(con, sql, params = list(chrom[i], start[i], end[i]))
      if (n > 1) r <- cbind(region_id = rep(sprintf("%s:%d-%d", chrom[i], start[i], end[i]), nrow(r)),
                            r, stringsAsFactors = FALSE)
      r
    })
    out <- do.call(rbind, res)
    rownames(out) <- NULL
    out
  }, db_file = db_file)
}

# ---- 6) Whole table / column names / build info --------------------------------
#' Whole table (about 3M rows; pass `columns` to keep memory down).
get_all_cpg_metadata <- function(columns = NULL, db_file = NULL, table_name = NULL) {
  table_name <- table_name %||% cpgdb_defaults()$table_name
  .with_con(function(con) {
    fields <- .check_table(con, table_name)
    cols <- if (is.null(columns)) fields else .resolve_columns(columns, fields)
    DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s",
                                 paste(DBI::dbQuoteIdentifier(con, cols), collapse = ", "),
                                 DBI::dbQuoteIdentifier(con, table_name)))
  }, db_file = db_file)
}

get_cpg_column_names <- function(prefix = NULL, include_cpg_id = FALSE,
                                 db_file = NULL, table_name = NULL) {
  table_name <- table_name %||% cpgdb_defaults()$table_name
  .with_con(function(con) {
    cols <- .check_table(con, table_name)
    if (!is.null(prefix)) {
      sel <- cols[startsWith(cols, prefix)]
      cols <- if (include_cpg_id && "cpg_id" %in% cols) c("cpg_id", sel) else sel
    }
    cols
  }, db_file = db_file)
}

get_anno_columns     <- function(include_cpg_id = FALSE, db_file = NULL, table_name = NULL)
  get_cpg_column_names("anno_", include_cpg_id, db_file, table_name)
get_var_columns      <- function(include_cpg_id = FALSE, db_file = NULL, table_name = NULL)
  get_cpg_column_names("var_", include_cpg_id, db_file, table_name)
get_arraymap_columns <- function(include_cpg_id = FALSE, db_file = NULL, table_name = NULL)
  get_cpg_column_names("arraymap_", include_cpg_id, db_file, table_name)

#' Build version, reference files and md5s, summary counts.
cpgdb_info <- function(db_file = NULL) {
  .with_con(function(con) {
    if (!DBI::dbExistsTable(con, "db_build_info")) stop("No db_build_info table (DB predates the pipeline).", call. = FALSE)
    DBI::dbReadTable(con, "db_build_info")
  }, db_file = db_file)
}


