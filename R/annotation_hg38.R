annotate_hg38_dmrs <- function(
    dmr_df,
    chr_col = "chr",
    start_col = "start_1based",
    end_col = "end_1based",
    dmr_id_col = NULL,
    include_intergenic_cpg = TRUE,
    extra_annotations = NULL
) {
  
  .require_optional(c("annotatr", "GenomicFeatures", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "AnnotationDbi"))
  required <- c(chr_col, start_col, end_col)
  
  if (!all(required %in% names(dmr_df))) {
    stop(
      "Missing required columns: ",
      paste(setdiff(required, names(dmr_df)), collapse = ", ")
    )
  }
  
  # ----------------------------------------------------------
  # Prepare DMR coordinates
  # ----------------------------------------------------------
  
  dmrs <- as.data.frame(dmr_df)
  
  dmrs$.dmr_index <- seq_len(nrow(dmrs))
  
  if (is.null(dmr_id_col)) {
    dmrs$.dmr_id <- paste0("DMR_", dmrs$.dmr_index)
  } else {
    if (!dmr_id_col %in% names(dmrs)) {
      stop("The dmr_id_col was not found in dmr_df.")
    }
    
    dmrs$.dmr_id <- as.character(dmrs[[dmr_id_col]])
  }
  
  dmrs$.chr <- as.character(dmrs[[chr_col]])
  dmrs$.chr <- ifelse(
    grepl("^chr", dmrs$.chr, ignore.case = TRUE),
    sub("^chr", "chr", dmrs$.chr, ignore.case = TRUE),
    paste0("chr", dmrs$.chr)
  )
  
  dmrs$.start <- as.integer(dmrs[[start_col]])
  dmrs$.end <- as.integer(dmrs[[end_col]])
  
  if (anyNA(dmrs[, c(".chr", ".start", ".end")])) {
    stop("Missing chromosome, start, or end values were found.")
  }
  
  if (any(dmrs$.start < 1L)) {
    stop(
      "Coordinates must be 1-based. ",
      "If using raw ENmix output, its 'start' column may be BED-like ",
      "and should be converted using start_1based = start + 1."
    )
  }
  
  if (any(dmrs$.start > dmrs$.end)) {
    stop("At least one DMR has start > end.")
  }
  
  gr_dmrs <- GenomicRanges::GRanges(
    seqnames = dmrs$.chr,
    ranges = IRanges::IRanges(
      start = dmrs$.start,
      end = dmrs$.end
    ),
    strand = "*",
    dmr_index = dmrs$.dmr_index,
    dmr_id = dmrs$.dmr_id
  )
  
  # Restrict to standard human chromosomes
  standard_chr <- paste0(
    "chr",
    c(1:22, "X", "Y", "M")
  )
  
  if (any(!as.character(GenomicRanges::seqnames(gr_dmrs)) %in% standard_chr)) {
    warning(
      "Some DMRs are not on standard chromosomes. ",
      "They will be retained but may not receive annotations."
    )
  }
  
  # ----------------------------------------------------------
  # Build annotatr annotations
  # ----------------------------------------------------------
  
  annotation_codes <- c(
    # Promoters, upstream regions, UTRs, exons and introns
    "hg38_basicgenes",
    
    # CpG islands, shores and shelves
    "hg38_cpgs",
    
    # FANTOM5 enhancers
    "hg38_enhancers_fantom",
    
    "hg38_genes_intergenic", "hg38_genes_intronexonboundaries"
    
  )
  
  annotations_gr <- annotatr::build_annotations(
    genome = "hg38",
    annotations = annotation_codes
  )
  
  # Optionally exclude the very broad inter-CpG-island category
  if (!include_intergenic_cpg) {
    annotations_gr <- annotations_gr[
      S4Vectors::mcols(annotations_gr)$type != "hg38_cpg_inter"
    ]
  }
  
  # Optional custom regulatory annotations, e.g. ENCODE cCREs
  if (!is.null(extra_annotations)) {
    if (!inherits(extra_annotations, "GRanges")) {
      stop("extra_annotations must be a GRanges object.")
    }
    
    if (is.null(S4Vectors::mcols(extra_annotations)$type)) {
      S4Vectors::mcols(extra_annotations)$type <- "hg38_custom_regulatory"
    }
    
    if (is.null(S4Vectors::mcols(extra_annotations)$id)) {
      S4Vectors::mcols(extra_annotations)$id <- paste0(
        "custom_",
        seq_along(extra_annotations)
      )
    }
    
    annotations_gr <- c(annotations_gr, extra_annotations)
  }
  
  # Harmonize chromosome naming
  GenomeInfoDb::seqlevelsStyle(annotations_gr) <- "UCSC"
  
  common_seqlevels <- intersect(
    GenomeInfoDb::seqlevels(gr_dmrs),
    GenomeInfoDb::seqlevels(annotations_gr)
  )
  
  gr_dmrs_overlap <- GenomeInfoDb::keepSeqlevels(
    gr_dmrs,
    common_seqlevels,
    pruning.mode = "coarse"
  )
  
  annotations_overlap <- GenomeInfoDb::keepSeqlevels(
    annotations_gr,
    common_seqlevels,
    pruning.mode = "coarse"
  )
  
  # ----------------------------------------------------------
  # Find overlapping annotations
  # ----------------------------------------------------------
  
  annotation_hits <- GenomicRanges::findOverlaps(
    gr_dmrs_overlap,
    annotations_overlap,
    ignore.strand = TRUE
  )
  
  if (length(annotation_hits) > 0) {
    
    annot_hit_df <- tibble::tibble(
      dmr_index = S4Vectors::mcols(gr_dmrs_overlap)$dmr_index[
        S4Vectors::queryHits(annotation_hits)
      ],
      annotation_type = as.character(
        S4Vectors::mcols(annotations_overlap)$type[
          S4Vectors::subjectHits(annotation_hits)
        ]
      ),
      annotation_id = as.character(
        S4Vectors::mcols(annotations_overlap)$id[
          S4Vectors::subjectHits(annotation_hits)
        ]
      ),
      gene_symbol = as.character(
        S4Vectors::mcols(annotations_overlap)$symbol[
          S4Vectors::subjectHits(annotation_hits)
        ]
      ),
      gene_id = as.character(
        S4Vectors::mcols(annotations_overlap)$gene_id[
          S4Vectors::subjectHits(annotation_hits)
        ]
      ),
      transcript_id = as.character(
        S4Vectors::mcols(annotations_overlap)$tx_id[
          S4Vectors::subjectHits(annotation_hits)
        ]
      )
    )
    
  } else {
    
    annot_hit_df <- tibble::tibble(
      dmr_index = integer(),
      annotation_type = character(),
      annotation_id = character(),
      gene_symbol = character(),
      gene_id = character(),
      transcript_id = character()
    )
  }
  
  collapse_unique <- function(x) {
    x <- unique(x[!is.na(x) & x != ""])
    
    if (length(x) == 0) {
      return(NA_character_)
    }
    
    paste(sort(x), collapse = ";")
  }
  
  # ----------------------------------------------------------
  # Collapse annotations to one row per DMR
  # ----------------------------------------------------------
  
  annotation_summary <- annot_hit_df %>%
    dplyr::group_by(dmr_index) %>%
    dplyr::summarise(
      overlapping_annotation_types = collapse_unique(
        sub("^hg38_", "", annotation_type)
      ),
      
      cpg_annotation = collapse_unique(
        sub(
          "^hg38_cpg_",
          "",
          annotation_type[grepl("^hg38_cpg_", annotation_type)]
        )
      ),
      
      gene_features = collapse_unique(
        sub(
          "^hg38_genes_",
          "",
          annotation_type[grepl("^hg38_genes_", annotation_type)]
        )
      ),
      
      overlapping_genes = collapse_unique(
        gene_symbol[
          grepl("^hg38_genes_", annotation_type)
        ]
      ),
      
      overlapping_gene_ids = collapse_unique(
        gene_id[
          grepl("^hg38_genes_", annotation_type)
        ]
      ),
      
      overlapping_transcripts = collapse_unique(
        transcript_id[
          grepl("^hg38_genes_", annotation_type)
        ]
      ),
      
      enhancer_overlap = any(
        grepl("enhancer", annotation_type, ignore.case = TRUE)
      ),
      
      enhancer_ids = collapse_unique(
        annotation_id[
          grepl("enhancer", annotation_type, ignore.case = TRUE)
        ]
      ),
      
      promoter_overlap = any(
        annotation_type == "hg38_genes_promoters"
      ),
      
      regulatory_annotations = collapse_unique(
        sub(
          "^hg38_",
          "",
          annotation_type[
            grepl(
              "promoter|enhancer|regulatory",
              annotation_type,
              ignore.case = TRUE
            )
          ]
        )
      ),
      
      .groups = "drop"
    )
  
  # ----------------------------------------------------------
  # Obtain gene-level hg38 knownGene loci
  # ----------------------------------------------------------
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  
  gene_gr <- GenomicFeatures::genes(
    txdb,
    single.strand.genes.only = TRUE
  )
  
  GenomeInfoDb::seqlevelsStyle(gene_gr) <- "UCSC"
  
  gene_gr <- GenomeInfoDb::keepStandardChromosomes(
    gene_gr,
    pruning.mode = "coarse"
  )
  
  gene_entrez <- names(gene_gr)
  
  gene_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_entrez,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  S4Vectors::mcols(gene_gr)$entrez_id <- gene_entrez
  S4Vectors::mcols(gene_gr)$symbol <- unname(gene_symbol)
  
  # ----------------------------------------------------------
  # Find closest gene
  # ----------------------------------------------------------
  
  closest_gene_index <- GenomicRanges::nearest(
    gr_dmrs,
    gene_gr,
    select = "arbitrary",
    ignore.strand = TRUE
  )
  
  closest_gene <- tibble::tibble(
    dmr_index = dmrs$.dmr_index,
    nearest_gene_symbol = NA_character_,
    nearest_gene_entrez = NA_character_,
    nearest_gene_chr = NA_character_,
    nearest_gene_start = NA_integer_,
    nearest_gene_end = NA_integer_,
    nearest_gene_strand = NA_character_,
    distance_to_gene_bp = NA_integer_,
    overlaps_nearest_gene_body = NA,
    nearest_gene_tss = NA_integer_,
    distance_to_tss_bp = NA_integer_,
    signed_distance_to_tss_bp = NA_integer_,
    position_relative_to_gene = NA_character_
  )
  
  valid <- which(!is.na(closest_gene_index))
  
  if (length(valid) > 0) {
    
    matched_genes <- gene_gr[closest_gene_index[valid]]
    matched_dmrs <- gr_dmrs[valid]
    
    # Edge-to-edge distance to gene body
    gene_body_distance <- pmax(
      BiocGenerics::start(matched_genes) - BiocGenerics::end(matched_dmrs),
      BiocGenerics::start(matched_dmrs) - BiocGenerics::end(matched_genes),
      0L
    )
    
    overlaps_gene <- IRanges::overlapsAny(
      matched_dmrs,
      matched_genes,
      ignore.strand = TRUE
    )
    
    # TSS is strand-specific
    gene_tss <- ifelse(
      as.character(GenomicRanges::strand(matched_genes)) == "-",
      BiocGenerics::end(matched_genes),
      BiocGenerics::start(matched_genes)
    )
    
    dmr_midpoint <- floor(
      (BiocGenerics::start(matched_dmrs) + BiocGenerics::end(matched_dmrs)) / 2
    )
    
    # Positive values are downstream relative to transcription
    signed_tss_distance <- ifelse(
      as.character(GenomicRanges::strand(matched_genes)) == "-",
      gene_tss - dmr_midpoint,
      dmr_midpoint - gene_tss
    )
    
    relative_position <- case_when(
      overlaps_gene ~ "overlaps_gene_body",
      signed_tss_distance < 0 ~ "upstream",
      signed_tss_distance > 0 ~ "downstream",
      TRUE ~ "at_TSS"
    )
    
    closest_gene[valid, ] <- tibble::tibble(
      dmr_index = dmrs$.dmr_index[valid],
      nearest_gene_symbol = as.character(
        S4Vectors::mcols(matched_genes)$symbol
      ),
      nearest_gene_entrez = as.character(
        S4Vectors::mcols(matched_genes)$entrez_id
      ),
      nearest_gene_chr = as.character(
        GenomicRanges::seqnames(matched_genes)
      ),
      nearest_gene_start = BiocGenerics::start(matched_genes),
      nearest_gene_end = BiocGenerics::end(matched_genes),
      nearest_gene_strand = as.character(
        GenomicRanges::strand(matched_genes)
      ),
      distance_to_gene_bp = as.integer(gene_body_distance),
      overlaps_nearest_gene_body = overlaps_gene,
      nearest_gene_tss = as.integer(gene_tss),
      distance_to_tss_bp = as.integer(
        abs(signed_tss_distance)
      ),
      signed_distance_to_tss_bp = as.integer(
        signed_tss_distance
      ),
      position_relative_to_gene = relative_position
    )
  }
  
  # ----------------------------------------------------------
  # Return one annotated row per DMR
  # ----------------------------------------------------------
  
  original_cols <- names(dmr_df)
  
  # Avoid returning the original ID column twice
  if (!is.null(dmr_id_col)) {
    original_cols <- setdiff(original_cols, dmr_id_col)
  }
  
  output <- dmrs %>%
    dplyr::left_join(
      annotation_summary,
      by = c(".dmr_index" = "dmr_index")
    ) %>%
    dplyr::left_join(
      closest_gene,
      by = c(".dmr_index" = "dmr_index")
    ) %>%
    dplyr::mutate(
      enhancer_overlap = dplyr::coalesce(
        enhancer_overlap,
        FALSE
      ),
      promoter_overlap = dplyr::coalesce(
        promoter_overlap,
        FALSE
      ),
      cpg_context = dplyr::case_when(
        grepl(
          "islands",
          cpg_annotation,
          ignore.case = TRUE
        ) ~ "CpG island",
        
        grepl(
          "shores",
          cpg_annotation,
          ignore.case = TRUE
        ) ~ "CpG shore",
        
        grepl(
          "shelves",
          cpg_annotation,
          ignore.case = TRUE
        ) ~ "CpG shelf",
        
        grepl(
          "inter",
          cpg_annotation,
          ignore.case = TRUE
        ) ~ "Inter-CGI",
        
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::rename(
      dmr_id = .dmr_id
    ) %>%
    dplyr::select(
      dmr_id,
      dplyr::all_of(original_cols),
      nearest_gene_symbol,
      nearest_gene_entrez,
      distance_to_gene_bp,
      overlaps_nearest_gene_body,
      position_relative_to_gene,
      nearest_gene_tss,
      distance_to_tss_bp,
      signed_distance_to_tss_bp,
      nearest_gene_chr,
      nearest_gene_start,
      nearest_gene_end,
      nearest_gene_strand,
      overlapping_genes,
      gene_features,
      promoter_overlap,
      enhancer_overlap,
      enhancer_ids,
      regulatory_annotations,
      cpg_context,
      cpg_annotation,
      overlapping_annotation_types,
      dplyr::everything()
    ) %>%
    dplyr::select(
      -dplyr::starts_with(".")
    )
  
  return(output)
}
