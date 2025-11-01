#' Generate summary reports for loci
#' 
#' @description
#' Create human-readable summary reports in markdown or text format
#' for genomic loci and their candidate genes.
#' 
#' @param locus_data Data frame with locus information
#' @param gene_data Data frame with gene annotations
#' @param output_format Format for output: "markdown", "text", "html"
#' @param include_sections Named list of sections to include
#' @param template Custom template for report generation
#' 
#' @return Character string with formatted report
#' 
#' @examples
#' locus <- data.frame(
#'   locus_id = "QTL_chr8",
#'   chr = 8,
#'   start = 28000000,
#'   end = 32000000,
#'   trait = "Heart weight"
#' )
#' genes <- data.frame(
#'   gene = c("Abcb10", "Acsl5"),
#'   biotype = c("protein_coding", "protein_coding")
#' )
#' report <- generate_locus_report(locus, genes)
#' cat(report)
#' 
#' @export
generate_locus_report <- function(locus_data,
                                 gene_data,
                                 output_format = "markdown",
                                 include_sections = NULL,
                                 template = NULL) {
  
  # Default sections
  if (is.null(include_sections)) {
    include_sections <- list(
      header = TRUE,
      summary_stats = TRUE,
      top_candidates = TRUE,
      gene_categories = TRUE,
      evidence_summary = TRUE,
      recommendations = TRUE
    )
  }
  
  # Initialize report
  report <- character()
  
  # Header section
  if (include_sections$header) {
    if (output_format == "markdown") {
      report <- c(report, 
                 sprintf("# Locus Report: %s", locus_data$locus_id[1]),
                 "",
                 sprintf("**Chromosome:** %s", locus_data$chr[1]),
                 sprintf("**Position:** %s - %s Mb", 
                        format(locus_data$start[1]/1e6, digits = 3),
                        format(locus_data$end[1]/1e6, digits = 3)),
                 sprintf("**Size:** %s Mb", 
                        format((locus_data$end[1] - locus_data$start[1])/1e6, digits = 3)),
                 "")
      
      if ("trait" %in% names(locus_data)) {
        report <- c(report,
                   sprintf("**Associated Trait:** %s", locus_data$trait[1]),
                   "")
      }
      
      report <- c(report, "---", "")
    } else {
      report <- c(report,
                 sprintf("LOCUS REPORT: %s", locus_data$locus_id[1]),
                 paste(rep("=", 50), collapse = ""),
                 sprintf("Chromosome: %s", locus_data$chr[1]),
                 sprintf("Position: %s - %s Mb", 
                        format(locus_data$start[1]/1e6, digits = 3),
                        format(locus_data$end[1]/1e6, digits = 3)),
                 "")
    }
  }
  
  # Summary statistics
  if (include_sections$summary_stats) {
    if (output_format == "markdown") {
      report <- c(report, "## Summary Statistics", "")
    } else {
      report <- c(report, "SUMMARY STATISTICS", paste(rep("-", 20), collapse = ""))
    }
    
    n_genes <- nrow(gene_data)
    n_protein <- sum(gene_data$biotype == "protein_coding", na.rm = TRUE)
    n_noncoding <- n_genes - n_protein
    
    stats_text <- c(
      sprintf("- Total genes: %d", n_genes),
      sprintf("- Protein-coding genes: %d (%.1f%%)", 
             n_protein, 100 * n_protein / n_genes),
      sprintf("- Non-coding genes: %d (%.1f%%)",
             n_noncoding, 100 * n_noncoding / n_genes)
    )
    
    # Add expression stats if available
    if (any(grepl("express|CPM|TPM", names(gene_data)))) {
      expr_col <- names(gene_data)[grepl("express|CPM|TPM", names(gene_data))][1]
      n_expressed <- sum(gene_data[[expr_col]] > 1, na.rm = TRUE)
      stats_text <- c(stats_text,
                     sprintf("- Expressed genes (>1 CPM): %d (%.1f%%)",
                            n_expressed, 100 * n_expressed / n_genes))
    }
    
    # Add orthology stats if available
    if ("human_ortholog" %in% names(gene_data)) {
      n_ortholog <- sum(!is.na(gene_data$human_ortholog))
      stats_text <- c(stats_text,
                     sprintf("- Genes with human orthologs: %d (%.1f%%)",
                            n_ortholog, 100 * n_ortholog / n_genes))
    }
    
    # Add variant stats if available
    if ("n_variants" %in% names(gene_data)) {
      n_with_var <- sum(gene_data$n_variants > 0, na.rm = TRUE)
      stats_text <- c(stats_text,
                     sprintf("- Genes with coding variants: %d (%.1f%%)",
                            n_with_var, 100 * n_with_var / n_genes))
    }
    
    # Add eQTL stats if available
    if ("has_eqtl" %in% names(gene_data)) {
      n_eqtl <- sum(gene_data$has_eqtl == TRUE, na.rm = TRUE)
      stats_text <- c(stats_text,
                     sprintf("- Genes with eQTL: %d (%.1f%%)",
                            n_eqtl, 100 * n_eqtl / n_genes))
    }
    
    report <- c(report, stats_text, "")
  }
  
  # Top candidates section
  if (include_sections$top_candidates) {
    if (output_format == "markdown") {
      report <- c(report, "## Top Candidate Genes", "")
    } else {
      report <- c(report, "TOP CANDIDATE GENES", paste(rep("-", 20), collapse = ""))
    }
    
    # Identify top candidates (simplified version)
    candidates <- identify_top_candidates(gene_data, n = 10)
    
    if (nrow(candidates) > 0) {
      if (output_format == "markdown") {
        report <- c(report, 
                   "| Gene | Biotype | Evidence | Priority Score |",
                   "|------|---------|----------|----------------|")
        
        for (i in 1:min(10, nrow(candidates))) {
          evidence <- character()
          if ("human_ortholog" %in% names(candidates) && !is.na(candidates$human_ortholog[i])) {
            evidence <- c(evidence, "Human ortholog")
          }
          if ("has_eqtl" %in% names(candidates) && candidates$has_eqtl[i]) {
            evidence <- c(evidence, "eQTL")
          }
          if ("n_variants" %in% names(candidates) && candidates$n_variants[i] > 0) {
            evidence <- c(evidence, sprintf("%d variants", candidates$n_variants[i]))
          }
          
          report <- c(report,
                     sprintf("| %s | %s | %s | %.1f |",
                            candidates$gene[i],
                            candidates$biotype[i],
                            paste(evidence, collapse = ", "),
                            candidates$priority_score[i]))
        }
        report <- c(report, "")
      } else {
        for (i in 1:min(10, nrow(candidates))) {
          report <- c(report,
                     sprintf("%d. %s (%s) - Score: %.1f",
                            i, 
                            candidates$gene[i],
                            candidates$biotype[i],
                            candidates$priority_score[i]))
        }
        report <- c(report, "")
      }
    }
  }
  
  # Gene categories section
  if (include_sections$gene_categories) {
    if (output_format == "markdown") {
      report <- c(report, "## Gene Categories", "")
    } else {
      report <- c(report, "GENE CATEGORIES", paste(rep("-", 20), collapse = ""))
    }
    
    # Categorize genes
    categories <- categorize_genes(gene_data)
    
    for (cat_name in names(categories)) {
      if (length(categories[[cat_name]]) > 0) {
        if (output_format == "markdown") {
          report <- c(report,
                     sprintf("### %s (%d genes)", cat_name, length(categories[[cat_name]])),
                     paste(categories[[cat_name]], collapse = ", "),
                     "")
        } else {
          report <- c(report,
                     sprintf("%s (%d genes):", cat_name, length(categories[[cat_name]])),
                     paste(" ", categories[[cat_name]], collapse = ", "),
                     "")
        }
      }
    }
  }
  
  # Evidence summary section
  if (include_sections$evidence_summary && "matched_keywords" %in% names(gene_data)) {
    if (output_format == "markdown") {
      report <- c(report, "## Phenotype Evidence", "")
    } else {
      report <- c(report, "PHENOTYPE EVIDENCE", paste(rep("-", 20), collapse = ""))
    }
    
    # Summarize phenotype matches
    genes_with_pheno <- gene_data[!is.na(gene_data$matched_keywords), ]
    if (nrow(genes_with_pheno) > 0) {
      # Count keyword frequencies
      all_keywords <- unlist(strsplit(genes_with_pheno$matched_keywords, "; "))
      keyword_freq <- table(all_keywords)
      keyword_freq <- sort(keyword_freq, decreasing = TRUE)
      
      if (output_format == "markdown") {
        report <- c(report, "**Most common phenotype associations:**", "")
        for (i in 1:min(5, length(keyword_freq))) {
          report <- c(report,
                     sprintf("- %s: %d genes", names(keyword_freq)[i], keyword_freq[i]))
        }
        report <- c(report, "")
      } else {
        report <- c(report, "Most common phenotype associations:")
        for (i in 1:min(5, length(keyword_freq))) {
          report <- c(report,
                     sprintf("  %s: %d genes", names(keyword_freq)[i], keyword_freq[i]))
        }
        report <- c(report, "")
      }
    }
  }
  
  # Recommendations section
  if (include_sections$recommendations) {
    if (output_format == "markdown") {
      report <- c(report, "## Recommendations for Follow-up", "")
    } else {
      report <- c(report, "RECOMMENDATIONS", paste(rep("-", 20), collapse = ""))
    }
    
    recommendations <- generate_recommendations(gene_data)
    
    if (output_format == "markdown") {
      report <- c(report, paste("1.", recommendations), "")
    } else {
      report <- c(report, recommendations, "")
    }
  }
  
  # Footer
  if (output_format == "markdown") {
    report <- c(report, 
               "---",
               sprintf("*Report generated on %s using PackRat*", Sys.Date()))
  } else {
    report <- c(report,
               paste(rep("=", 50), collapse = ""),
               sprintf("Report generated on %s using PackRat", Sys.Date()))
  }
  
  # Combine report lines
  if (output_format == "html") {
    # Convert markdown to HTML (basic conversion)
    html_report <- paste(report, collapse = "\n")
    html_report <- gsub("^# (.+)$", "<h1>\\1</h1>", html_report)
    html_report <- gsub("^## (.+)$", "<h2>\\1</h2>", html_report)
    html_report <- gsub("^### (.+)$", "<h3>\\1</h3>", html_report)
    html_report <- gsub("\\*\\*(.+?)\\*\\*", "<strong>\\1</strong>", html_report)
    html_report <- gsub("^- (.+)$", "<li>\\1</li>", html_report)
    return(paste("<html><body>", html_report, "</body></html>"))
  } else {
    return(paste(report, collapse = "\n"))
  }
}

#' Identify top candidate genes
#' 
#' @description
#' Internal function to identify and score top candidate genes.
#' 
#' @param gene_data Gene data frame
#' @param n Number of top candidates to return
#' @return Data frame of top candidates with scores
#' 
#' @keywords internal
identify_top_candidates <- function(gene_data, n = 10) {
  
  # Copy data to avoid modifying original
  candidates <- gene_data
  
  # Initialize priority score
  candidates$priority_score <- 0
  
  # Scoring based on available evidence
  if ("biotype" %in% names(candidates)) {
    candidates$priority_score <- candidates$priority_score + 
      ifelse(candidates$biotype == "protein_coding", 3, 0)
  }
  
  if ("human_ortholog" %in% names(candidates)) {
    candidates$priority_score <- candidates$priority_score + 
      ifelse(!is.na(candidates$human_ortholog), 2, 0)
  }
  
  if (any(grepl("express|CPM|TPM", names(candidates)))) {
    expr_col <- names(candidates)[grepl("express|CPM|TPM", names(candidates))][1]
    expr_median <- median(candidates[[expr_col]], na.rm = TRUE)
    candidates$priority_score <- candidates$priority_score + 
      ifelse(candidates[[expr_col]] > expr_median, 1, 0) +
      ifelse(candidates[[expr_col]] > 2 * expr_median, 1, 0)
  }
  
  if ("has_eqtl" %in% names(candidates)) {
    candidates$priority_score <- candidates$priority_score + 
      ifelse(candidates$has_eqtl == TRUE, 2, 0)
  }
  
  if ("n_variants" %in% names(candidates)) {
    candidates$priority_score <- candidates$priority_score + 
      ifelse(candidates$n_variants > 0, 1, 0) +
      ifelse(candidates$n_variants > 2, 1, 0)
  }
  
  if ("n_phenotypes" %in% names(candidates)) {
    candidates$priority_score <- candidates$priority_score + 
      ifelse(candidates$n_phenotypes > 0, 1, 0) +
      ifelse(candidates$n_phenotypes > 3, 1, 0)
  }
  
  # Sort by score and return top n
  candidates <- candidates[order(candidates$priority_score, decreasing = TRUE), ]
  return(head(candidates, n))
}

#' Categorize genes by biotype and evidence
#' 
#' @description
#' Internal function to categorize genes for reporting.
#' 
#' @param gene_data Gene data frame
#' @return Named list of gene categories
#'
#' @keywords internal
#' @importFrom stats quantile
categorize_genes <- function(gene_data) {
  
  categories <- list()
  
  # Biotype categories
  if ("biotype" %in% names(gene_data)) {
    categories[["Protein-coding"]] <- gene_data$gene[gene_data$biotype == "protein_coding"]
    categories[["lncRNA"]] <- gene_data$gene[grepl("lncRNA|lincRNA", gene_data$biotype)]
    categories[["miRNA"]] <- gene_data$gene[gene_data$biotype == "miRNA"]
  }
  
  # Evidence-based categories
  if ("has_eqtl" %in% names(gene_data)) {
    categories[["With eQTL"]] <- gene_data$gene[gene_data$has_eqtl == TRUE]
  }
  
  if ("n_variants" %in% names(gene_data)) {
    categories[["With coding variants"]] <- gene_data$gene[gene_data$n_variants > 0]
  }
  
  # Expression categories
  if (any(grepl("express|CPM|TPM", names(gene_data)))) {
    expr_col <- names(gene_data)[grepl("express|CPM|TPM", names(gene_data))][1]
    expr_q75 <- quantile(gene_data[[expr_col]], 0.75, na.rm = TRUE)
    categories[["Highly expressed"]] <- gene_data$gene[gene_data[[expr_col]] > expr_q75]
  }
  
  return(categories)
}

#' Generate recommendations based on locus analysis
#' 
#' @description
#' Internal function to generate follow-up recommendations.
#' 
#' @param gene_data Gene data frame
#' @return Character vector of recommendations
#' 
#' @keywords internal
generate_recommendations <- function(gene_data) {
  
  recommendations <- character()
  
  # Check data completeness
  n_genes <- nrow(gene_data)
  n_protein <- sum(gene_data$biotype == "protein_coding", na.rm = TRUE)
  
  # High-priority recommendations
  if ("has_eqtl" %in% names(gene_data)) {
    n_eqtl <- sum(gene_data$has_eqtl == TRUE, na.rm = TRUE)
    if (n_eqtl > 0 && n_eqtl <= 5) {
      recommendations <- c(recommendations,
                          sprintf("Focus on %d genes with eQTL evidence for functional validation", n_eqtl))
    }
  }
  
  if ("n_variants" %in% names(gene_data)) {
    genes_with_nonsense <- gene_data$gene[grepl("nonsense|stop", gene_data$variant_type, ignore.case = TRUE)]
    if (length(genes_with_nonsense) > 0 && length(genes_with_nonsense) <= 3) {
      recommendations <- c(recommendations,
                          sprintf("Prioritize genes with nonsense variants: %s", 
                                 paste(genes_with_nonsense, collapse = ", ")))
    }
  }
  
  # General recommendations
  if (n_protein > 20) {
    recommendations <- c(recommendations,
                        "Consider additional filtering criteria to narrow candidate list")
  }
  
  if (!any(grepl("express|CPM|TPM", names(gene_data)))) {
    recommendations <- c(recommendations,
                        "Add tissue-specific expression data to improve prioritization")
  }
  
  if (!"human_ortholog" %in% names(gene_data)) {
    recommendations <- c(recommendations,
                        "Include human orthology information for translational relevance")
  }
  
  if (length(recommendations) == 0) {
    recommendations <- "Review top candidates and consider functional validation studies"
  }
  
  return(recommendations)
}