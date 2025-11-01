# Declare global variables for data.table NSE
utils::globalVariables(c(".", ":=", ".N", "chr", "pos", "p_value",
                         "variant_type", "n_variants", "human_gene_id",
                         "biotype", "locus_id", "n_genes"))

#' Join gene tables with optional annotations
#'
#' @description
#' Core function for joining gene tables with various annotation sources.
#' Follows relational database principles to combine gene-level data from
#' multiple sources into a unified table.
#' 
#' @param genes_df Data frame with gene identifiers (required)
#' @param gene_id_col Name of gene ID column in genes_df (default: "gene_id")
#' @param gene_symbol_col Name of gene symbol column (default: "gene_symbol")
#' @param orthology_df Optional orthology mapping table
#' @param expression_df Optional expression data table
#' @param eqtl_df Optional eQTL results table
#' @param variants_df Optional variant annotation table
#' @param phenotypes_df Optional phenotype associations table
#' @param join_type Type of join to use ("left", "inner", "full")
#' 
#' @return A data.table with joined gene information
#' 
#' @examples
#' # Basic usage with just genes
#' genes <- data.frame(gene_id = c("ENSMUSG001", "ENSMUSG002"),
#'                     gene_symbol = c("Abcb10", "Acsl5"))
#' result <- build_gene_table(genes)
#' 
#' # Add expression data
#' expr <- data.frame(gene_id = c("ENSMUSG001", "ENSMUSG002"),
#'                    CPM_heart = c(150.2, 89.5))
#' result <- build_gene_table(genes, expression_df = expr)
#' 
#' @export
#' @importFrom data.table data.table setDT merge.data.table := .N copy
#' @importFrom dplyr left_join inner_join full_join

build_gene_table <- function(genes_df,
                           gene_id_col = "gene_id",
                           gene_symbol_col = "gene_symbol",
                           orthology_df = NULL,
                           expression_df = NULL,
                           eqtl_df = NULL,
                           variants_df = NULL,
                           phenotypes_df = NULL,
                           join_type = "left") {
  
  # Convert to data.table for efficiency
  genes_dt <- data.table::as.data.table(genes_df)
  
  # Standardize column names
  if (gene_id_col != "gene_id") {
    data.table::setnames(genes_dt, gene_id_col, "gene_id")
  }
  if (gene_symbol_col != "gene_symbol" && gene_symbol_col %in% names(genes_dt)) {
    data.table::setnames(genes_dt, gene_symbol_col, "gene_symbol")
  }
  
  # Join orthology if provided
  if (!is.null(orthology_df)) {
    orth_dt <- data.table::as.data.table(orthology_df)
    genes_dt <- merge(genes_dt, orth_dt, 
                     by = "gene_id", 
                     all.x = (join_type != "inner"),
                     all.y = (join_type == "full"))
  }
  
  # Join expression data if provided
  if (!is.null(expression_df)) {
    expr_dt <- data.table::as.data.table(expression_df)
    genes_dt <- merge(genes_dt, expr_dt, 
                     by = "gene_id", 
                     all.x = (join_type != "inner"),
                     all.y = (join_type == "full"))
  }
  
  # Join eQTL data if provided
  if (!is.null(eqtl_df)) {
    eqtl_dt <- data.table::as.data.table(eqtl_df)
    # Summarize eQTL by gene if multiple per gene
    if ("p_value" %in% names(eqtl_dt)) {
      eqtl_summary <- eqtl_dt[, .(min_eqtl_p = min(p_value),
                                  n_eqtl = .N,
                                  has_eqtl = TRUE), 
                             by = gene_id]
      genes_dt <- merge(genes_dt, eqtl_summary, 
                       by = "gene_id", 
                       all.x = TRUE)
      genes_dt[is.na(has_eqtl), has_eqtl := FALSE]
    } else {
      genes_dt <- merge(genes_dt, eqtl_dt, 
                       by = "gene_id", 
                       all.x = (join_type != "inner"),
                       all.y = (join_type == "full"))
    }
  }
  
  # Join variant data if provided
  if (!is.null(variants_df)) {
    var_dt <- data.table::as.data.table(variants_df)
    # Summarize variants by gene if multiple per gene
    if ("variant_type" %in% names(var_dt)) {
      var_summary <- var_dt[, .(n_variants = .N,
                               has_missense = any(variant_type == "missense"),
                               has_nonsense = any(variant_type == "nonsense")), 
                           by = gene_id]
      genes_dt <- merge(genes_dt, var_summary, 
                       by = "gene_id", 
                       all.x = TRUE)
      genes_dt[is.na(n_variants), n_variants := 0]
    } else {
      genes_dt <- merge(genes_dt, var_dt, 
                       by = "gene_id", 
                       all.x = (join_type != "inner"),
                       all.y = (join_type == "full"))
    }
  }
  
  # Join phenotype data if provided
  if (!is.null(phenotypes_df)) {
    pheno_dt <- data.table::as.data.table(phenotypes_df)
    genes_dt <- merge(genes_dt, pheno_dt, 
                     by = "gene_id", 
                     all.x = (join_type != "inner"),
                     all.y = (join_type == "full"))
  }
  
  return(genes_dt)
}

#' Aggregate gene data by locus
#' 
#' @description
#' Summarizes gene-level data within genomic loci, useful for QTL/GWAS regions.
#' 
#' @param genes_dt Gene table from build_gene_table()
#' @param locus_ranges GRanges or data frame with chr, start, end columns
#' @param gene_chr_col Column name for chromosome in genes_dt
#' @param gene_start_col Column name for gene start position
#' @param gene_end_col Column name for gene end position
#' 
#' @return Data table with genes assigned to loci
#'
#' @export
#' @importFrom data.table data.table setDT := .N
aggregate_genes_by_locus <- function(genes_dt,
                                    locus_ranges,
                                    gene_chr_col = "chromosome",
                                    gene_start_col = "start",
                                    gene_end_col = "end") {
  
  # Convert locus ranges to data.table if needed
  if (inherits(locus_ranges, "GRanges")) {
    loci_dt <- data.table::data.table(
      locus_id = names(locus_ranges),
      chr = as.character(GenomeInfoDb::seqnames(locus_ranges)),
      start = BiocGenerics::start(locus_ranges),
      end = BiocGenerics::end(locus_ranges)
    )
  } else {
    loci_dt <- data.table::as.data.table(locus_ranges)
  }
  
  # Standardize column names
  if (gene_chr_col != "chr") {
    data.table::setnames(genes_dt, gene_chr_col, "chr", skip_absent = TRUE)
  }
  if (gene_start_col != "start") {
    data.table::setnames(genes_dt, gene_start_col, "start", skip_absent = TRUE)
  }
  if (gene_end_col != "end") {
    data.table::setnames(genes_dt, gene_end_col, "end", skip_absent = TRUE)
  }
  
  # Non-equi join to find overlapping genes
  genes_dt[, gene_mid := (start + end) / 2]
  
  # For each locus, find overlapping genes
  result <- loci_dt[genes_dt, 
                   on = .(chr = chr, start <= gene_mid, end >= gene_mid),
                   nomatch = 0,
                   allow.cartesian = TRUE]
  
  # Add locus summary statistics
  locus_summary <- result[, .(
    n_genes = .N,
    n_protein_coding = sum(biotype == "protein_coding", na.rm = TRUE),
    n_with_expression = sum(!is.na(CPM_heart), na.rm = TRUE),
    n_with_eqtl = sum(has_eqtl == TRUE, na.rm = TRUE),
    n_with_variants = sum(n_variants > 0, na.rm = TRUE)
  ), by = locus_id]
  
  result <- merge(result, locus_summary, by = "locus_id")
  
  return(result)
}

#' Filter genes by criteria
#' 
#' @description
#' Apply standardized filtering criteria to gene tables.
#' 
#' @param genes_dt Gene table from build_gene_table()
#' @param require_human_ortholog Require human ortholog (default: FALSE)
#' @param require_protein_coding Require protein coding genes (default: FALSE)
#' @param min_expression Minimum expression threshold (default: NULL)
#' @param require_eqtl Require eQTL association (default: FALSE)
#' @param require_variants Require coding variants (default: FALSE)
#' 
#' @return Filtered gene data table
#'
#' @export
#' @importFrom data.table copy
filter_genes <- function(genes_dt,
                        require_human_ortholog = FALSE,
                        require_protein_coding = FALSE,
                        min_expression = NULL,
                        require_eqtl = FALSE,
                        require_variants = FALSE) {
  
  result <- copy(genes_dt)
  
  if (require_human_ortholog && "human_gene_id" %in% names(result)) {
    result <- result[!is.na(human_gene_id)]
  }
  
  if (require_protein_coding && "biotype" %in% names(result)) {
    result <- result[biotype == "protein_coding"]
  }
  
  if (!is.null(min_expression) && "CPM_heart" %in% names(result)) {
    result <- result[CPM_heart >= min_expression]
  }
  
  if (require_eqtl && "has_eqtl" %in% names(result)) {
    result <- result[has_eqtl == TRUE]
  }
  
  if (require_variants && "n_variants" %in% names(result)) {
    result <- result[n_variants > 0]
  }
  
  return(result)
}