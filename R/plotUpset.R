#' Generate an UpSet Plot from locusPackRat Project Data
#'
#' Creates an UpSet plot visualization showing set intersections for genes
#' based on criteria derived from supplementary tables in a locusPackRat project.
#'
#' @param criteria Named list of criterion specifications. Each element should be
#'   a list with: \code{table} (table name), \code{column} (column to evaluate),
#'   and \code{condition} (expression string). See Details.
#' @param project_dir Character: Path to the locusPackRat project directory
#'   containing the \code{.locusPackRat} folder.
#' @param highlight_genes Character vector (optional): Gene symbols to highlight
#'   with a distinct color in the plot.
#' @param highlight_intersections List (optional): List of character vectors specifying
#'   which intersections to highlight. Each vector names the criteria that must
#'   all be TRUE for that intersection.
#' @param highlight_color Character: Color for highlighted genes (default "#a03b60").
#' @param intersection_color Character: Color for highlighted intersections (default "#53add2").
#' @param min_size Integer: Minimum intersection size to display (default 1).
#' @param sort_intersections Character: How to sort intersections - "descending",
#'   "ascending", or "cardinality" (default "descending").
#' @param sort_sets Character: How to sort sets - "descending", "ascending" (default "descending").
#' @param width Numeric: Width of output in inches (default 8).
#' @param height Numeric: Height of output in inches (default 5).
#' @param text_scale Numeric (optional): Scaling factor for text elements.
#'   If \code{NULL}, a default scaling based on width/height is applied.
#' @param font_family Character (optional): Font family for all text.
#' @param output_file Character: Output file name (default "upset_plot.pdf").
#' @param device Character: Output device - "pdf", "png", or "svg" (default "pdf").
#' @param dpi Integer: Resolution for PNG output (default 300).
#'
#' @return Invisible ggplot object (plot is also saved to disk).
#'
#' @details
#' ## Criteria Specification
#'
#' The \code{criteria} parameter accepts a named list where each element defines one
#' criterion (set) in the UpSet plot. Each criterion is a list with:
#'
#' \itemize{
#'   \item \code{table}: Name of supplementary table (without .csv)
#'   \item \code{column}: Column name to evaluate
#'   \item \code{condition}: One of:
#'     \itemize{
#'       \item A comparison operator with value: \code{"< 0.05"}, \code{"== TRUE"}, \code{"!= NA"}
#'       \item \code{"exists"} - TRUE if gene has any row in this table
#'       \item \code{"not_na"} - TRUE if column value is not NA
#'       \item \code{"matches <pattern>"} - TRUE if value matches regex pattern
#'     }
#' }
#'
#' @examples
#' \dontrun{
#' # Define criteria
#' criteria <- list(
#'   "DE (FDR < 0.05)" = list(
#'     table = "deseq_results",
#'     column = "padj",
#'     condition = "< 0.05"
#'   ),
#'   "Has cis-eQTL" = list(
#'     table = "eqtl_data",
#'     column = "gene_symbol",
#'     condition = "exists"
#'   )
#' )
#'
#' # Generate plot
#' generateUpsetPlot(
#'   criteria = criteria,
#'   project_dir = "my_analysis",
#'   highlight_genes = c("Tspan5", "Pdlim5"),
#'   output_file = "upset_plot.pdf"
#' )
#' }
#'
#' @importFrom data.table fread setDT as.data.table
#' @importFrom jsonlite read_json
#' @importFrom ComplexUpset upset intersection_size upset_modify_themes upset_query
#' @importFrom ggplot2 ggsave theme element_text element_blank labs
#'
#' @export
generateUpsetPlot <- function(
    criteria,
    project_dir = ".",
    highlight_genes = NULL,
    highlight_intersections = NULL,
    highlight_color = "#a03b60",
    intersection_color = "#53add2",
    min_size = 1,
    sort_intersections = c("descending", "ascending", "cardinality"),
    sort_sets = c("descending", "ascending"),
    width = 8,
    height = 5,
    text_scale = NULL,
    font_family = NULL,
    output_file = "upset_plot.pdf",
    device = c("pdf", "png", "svg"),
    dpi = 300
) {
  # Input validation
  sort_intersections <- match.arg(sort_intersections)
  sort_sets <- match.arg(sort_sets)
  device <- match.arg(device)

  if (!is.list(criteria) || length(criteria) == 0) {
    stop("criteria must be a non-empty named list")
  }
  if (is.null(names(criteria)) || any(names(criteria) == "")) {
    stop("All criteria must be named")
  }

  # Project validation
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    stop("No .locusPackRat directory found. Run initPackRat() first.")
  }

  config <- jsonlite::read_json(file.path(packrat_dir, "config.json"))

  # Text scaling (follows plotZoom pattern)
  if (is.null(text_scale)) {
    text_scale <- min(width / 8, height / 5)
  }

  # Get gene universe
  gene_universe <- .getGeneUniverse(packrat_dir, config)
  message(sprintf("Gene universe: %d genes", length(gene_universe)))

  # Build boolean matrix from criteria
  message("Building upset matrix from criteria...")
  upset_data <- .buildUpsetMatrix(
    criteria = criteria,
    gene_universe = gene_universe,
    packrat_dir = packrat_dir
  )

  message(sprintf("Built upset matrix: %d genes x %d criteria",
                  nrow(upset_data), length(criteria)))

  # Build highlight queries
  queries <- .buildUpsetQueries(
    highlight_genes = highlight_genes,
    highlight_intersections = highlight_intersections,
    highlight_color = highlight_color,
    intersection_color = intersection_color,
    upset_data = upset_data
  )

  # Build intersection panel
  inter_panel <- ComplexUpset::intersection_size(
    counts = TRUE,
    text = list(vjust = -0.6, size = 3 * text_scale)
  ) + ggplot2::labs(y = "Intersection size")

  # Build themes
  theme_mods <- list(
    'intersections_matrix' = ggplot2::theme(
      axis.title.y = ggplot2::element_blank()
    ),
    'overall_sizes' = ggplot2::theme(
      axis.text = ggplot2::element_text(angle = 90)
    )
  )

  # Add font family if specified
  if (!is.null(font_family)) {
    theme_mods[['default']] <- ggplot2::theme(
      text = ggplot2::element_text(family = font_family)
    )
  }

  upset_themes <- ComplexUpset::upset_modify_themes(theme_mods)

  # Create plot
  criteria_names <- names(criteria)

  p <- ComplexUpset::upset(
    upset_data,
    intersect = criteria_names,
    name = "Criteria",
    sort_intersections = sort_intersections,
    sort_sets = sort_sets,
    width_ratio = 0.25,
    min_size = min_size,
    base_annotations = list("Intersection size" = inter_panel),
    themes = upset_themes,
    queries = queries
  )

  # Ensure output directory exists
  output_dir <- file.path(packrat_dir, "output")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save output
  output_path <- file.path(output_dir, output_file)

  if (device == "pdf") {
    ggplot2::ggsave(
      output_path,
      p,
      width = width,
      height = height,
      device = grDevices::cairo_pdf
    )
  } else {
    ggplot2::ggsave(
      output_path,
      p,
      width = width,
      height = height,
      device = device,
      dpi = dpi
    )
  }

  message(sprintf("Saved upset plot to %s", output_path))
  invisible(p)
}


#' Check Criteria Validity for UpSet Plot
#'
#' Validates that criteria specifications reference existing tables/columns
#' and that conditions parse correctly.
#'
#' @param criteria Named list of criterion specifications (same format as generateUpsetPlot).
#' @param project_dir Character: Path to project directory.
#' @param verbose Logical: Print detailed validation results (default TRUE).
#'
#' @return Invisible TRUE if all criteria valid. Prints validation results.
#'
#' @examples
#' \dontrun{
#' criteria <- list(
#'   "DE Genes" = list(table = "deseq", column = "padj", condition = "< 0.05")
#' )
#' checkCriteria(criteria, project_dir = "my_analysis")
#' }
#'
#' @importFrom data.table fread
#'
#' @export
checkCriteria <- function(criteria, project_dir = ".", verbose = TRUE) {
  # Validate inputs
  if (!is.list(criteria) || length(criteria) == 0) {
    stop("criteria must be a non-empty named list")
  }
  if (is.null(names(criteria)) || any(names(criteria) == "")) {
    stop("All criteria must be named")
  }

  # Project validation
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    stop("No .locusPackRat directory found. Run initPackRat() first.")
  }

  config <- jsonlite::read_json(file.path(packrat_dir, "config.json"))
  supp_dir <- file.path(packrat_dir, "supplementary")

  # Get gene universe for count checking
  gene_universe <- .getGeneUniverse(packrat_dir, config)

  if (verbose) {
    message(sprintf("Checking %d criteria against project tables...", length(criteria)))
    message(sprintf("Gene universe: %d genes\n", length(gene_universe)))
  }

  n_passed <- 0
  n_failed <- 0

  for (crit_name in names(criteria)) {
    crit <- criteria[[crit_name]]

    # Check structure
    if (!all(c("table", "column", "condition") %in% names(crit))) {
      if (verbose) {
        message(sprintf("  x \"%s\": Missing required fields (need: table, column, condition)",
                        crit_name))
      }
      n_failed <- n_failed + 1
      next
    }

    # Check table exists
    table_path <- file.path(supp_dir, paste0(crit$table, ".csv"))
    if (!file.exists(table_path)) {
      if (verbose) {
        message(sprintf("  x \"%s\": table '%s' NOT FOUND", crit_name, crit$table))
      }
      n_failed <- n_failed + 1
      next
    }

    # Load table header to check column
    table_header <- data.table::fread(table_path, nrows = 0)

    # For "exists" condition, column can be the gene column itself
    if (crit$condition != "exists" && !crit$column %in% names(table_header)) {
      if (verbose) {
        message(sprintf("  x \"%s\": column '%s' not found in table '%s'",
                        crit_name, crit$column, crit$table))
        message(sprintf("      Available columns: %s",
                        paste(names(table_header), collapse = ", ")))
      }
      n_failed <- n_failed + 1
      next
    }

    # Load full table for evaluation
    table_data <- data.table::fread(table_path)

    # Get column type
    if (crit$column %in% names(table_data)) {
      col_type <- class(table_data[[crit$column]])[1]
    } else {
      col_type <- "N/A"
    }

    # Try to evaluate condition
    tryCatch({
      bool_vec <- .evaluateCriterion(
        table_data = table_data,
        column = crit$column,
        condition = crit$condition,
        gene_universe = gene_universe
      )

      n_match <- sum(bool_vec)
      pct_match <- round(100 * mean(bool_vec), 1)

      if (verbose) {
        message(sprintf("  v \"%s\": table '%s' found, column '%s' exists (%s)",
                        crit_name, crit$table, crit$column, col_type))
        message(sprintf("      Condition '%s' valid - %d genes match (%.1f%%)",
                        crit$condition, n_match, pct_match))
      }
      n_passed <- n_passed + 1

    }, error = function(e) {
      if (verbose) {
        message(sprintf("  x \"%s\": condition '%s' failed to evaluate",
                        crit_name, crit$condition))
        message(sprintf("      Error: %s", e$message))
      }
      n_failed <<- n_failed + 1
    })
  }

  if (verbose) {
    message("")
    if (n_failed == 0) {
      message(sprintf("All %d criteria validated successfully!", n_passed))
    } else {
      message(sprintf("Validation complete: %d passed, %d failed", n_passed, n_failed))
    }
  }

  if (n_failed > 0) {
    stop(sprintf("%d criteria failed validation", n_failed))
  }

  invisible(TRUE)
}


#' Get gene universe from project
#' @noRd
.getGeneUniverse <- function(packrat_dir, config) {
  if (config$mode == "gene") {
    genes_path <- file.path(packrat_dir, "input", "genes.csv")
    if (!file.exists(genes_path)) {
      stop("genes.csv not found in project input directory")
    }
    genes_dt <- data.table::fread(genes_path)

    # Find gene column
    gene_col <- intersect(
      c("gene_symbol", "gene", "mouse_gene_symbol", "human_gene_symbol"),
      names(genes_dt)
    )[1]

    if (is.na(gene_col)) {
      stop("Could not find gene identifier column in genes.csv")
    }

    return(unique(genes_dt[[gene_col]]))

  } else {
    # Region mode - extract genes from regions
    regions_path <- file.path(packrat_dir, "input", "regions.csv")
    if (!file.exists(regions_path)) {
      stop("regions.csv not found in project input directory")
    }
    regions_dt <- data.table::fread(regions_path)

    if (!"genes" %in% names(regions_dt)) {
      stop("regions.csv must have a 'genes' column")
    }

    # Parse comma-separated genes
    all_genes <- unlist(strsplit(regions_dt$genes, ",\\s*"))
    return(unique(all_genes[all_genes != "" & !is.na(all_genes)]))
  }
}


#' Build boolean matrix from criteria specifications
#' @noRd
.buildUpsetMatrix <- function(criteria, gene_universe, packrat_dir) {
  # Start with gene universe
  upset_dt <- data.table::data.table(gene = gene_universe)

  supp_dir <- file.path(packrat_dir, "supplementary")

  for (crit_name in names(criteria)) {
    crit <- criteria[[crit_name]]

    # Load table
    table_path <- file.path(supp_dir, paste0(crit$table, ".csv"))
    table_data <- data.table::fread(table_path)

    # Evaluate criterion
    bool_vec <- .evaluateCriterion(
      table_data = table_data,
      column = crit$column,
      condition = crit$condition,
      gene_universe = gene_universe
    )

    # Add to matrix
    upset_dt[, (crit_name) := bool_vec]

    message(sprintf("  %s: %d genes TRUE (%.1f%%)",
                    crit_name, sum(bool_vec), 100 * mean(bool_vec)))
  }

  return(as.data.frame(upset_dt))
}


#' Evaluate a single criterion to produce boolean vector
#' @noRd
.evaluateCriterion <- function(table_data, column, condition, gene_universe) {
  # Detect gene column in table
  gene_col <- intersect(
    c("gene_symbol", "gene", "mouse_gene_symbol", "human_gene_symbol"),
    names(table_data)
  )[1]

  if (is.na(gene_col)) {
    stop("Could not find gene identifier column in table")
  }

  # Parse condition
  condition <- trimws(condition)

  if (condition == "exists") {
    # Gene exists in table
    return(gene_universe %in% table_data[[gene_col]])

  } else if (condition == "not_na") {
    # Column value is not NA for this gene
    genes_with_value <- table_data[!is.na(get(column)), get(gene_col)]
    return(gene_universe %in% genes_with_value)

  } else if (grepl("^matches ", condition, ignore.case = TRUE)) {
    # Regex match
    pattern <- sub("^matches\\s+", "", condition, ignore.case = TRUE)
    matching_rows <- table_data[grepl(pattern, get(column), ignore.case = TRUE)]
    matching_genes <- matching_rows[[gene_col]]
    return(gene_universe %in% matching_genes)

  } else {
    # Comparison operator: <, >, <=, >=, ==, !=, %in%
    # Build expression and evaluate
    expr_str <- sprintf("%s %s", column, condition)

    tryCatch({
      matching_rows <- table_data[eval(parse(text = expr_str))]
      matching_genes <- matching_rows[[gene_col]]
      return(gene_universe %in% matching_genes)
    }, error = function(e) {
      stop(sprintf("Failed to evaluate condition '%s': %s", condition, e$message))
    })
  }
}


#' Build ComplexUpset queries for highlighting
#' @noRd
.buildUpsetQueries <- function(highlight_genes, highlight_intersections,
                                highlight_color, intersection_color, upset_data) {
  queries <- list()

  # Gene highlighting - convert to intersection highlighting
  if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
    # Filter to genes that exist in the data
    valid_genes <- intersect(highlight_genes, upset_data$gene)

    if (length(valid_genes) > 0) {
      # Get criteria columns (all columns except 'gene')
      criteria_cols <- setdiff(names(upset_data), "gene")

      # For each highlighted gene, find its intersection signature
      for (g in valid_genes) {
        gene_row <- upset_data[upset_data$gene == g, , drop = FALSE]
        # Get which criteria are TRUE for this gene
        true_criteria <- criteria_cols[unlist(gene_row[, criteria_cols])]

        if (length(true_criteria) > 0) {
          queries[[length(queries) + 1]] <- ComplexUpset::upset_query(
            intersect = true_criteria,
            color = highlight_color,
            fill = highlight_color,
            only_components = c('intersections_matrix', 'Intersection size')
          )
        }
      }
      message(sprintf("Highlighting %d genes via their intersections", length(valid_genes)))
    } else {
      warning("None of the specified highlight_genes found in data")
    }
  }

  # Intersection highlighting
  if (!is.null(highlight_intersections)) {
    for (i in seq_along(highlight_intersections)) {
      intersect_spec <- highlight_intersections[[i]]

      # Validate intersection criteria exist
      missing <- setdiff(intersect_spec, names(upset_data))
      if (length(missing) > 0) {
        warning(sprintf("Intersection %d: criteria not found: %s",
                        i, paste(missing, collapse = ", ")))
        next
      }

      queries[[length(queries) + 1]] <- ComplexUpset::upset_query(
        intersect = intersect_spec,
        color = intersection_color,
        fill = intersection_color,
        only_components = c('intersections_matrix')
      )
    }
    message(sprintf("Highlighting %d intersections", length(highlight_intersections)))
  }

  if (length(queries) == 0) return(NULL)
  return(queries)
}
