#' Initialize locusPackRat Project
#'
#' Creates a .locusPackRat directory structure and initializes project with gene or region data
#'
#' @param data data.frame or data.table with gene/region information
#' @param mode Character: "gene" or "region" mode
#' @param species Character: "human" or "mouse"
#' @param genome Character: "hg19", "hg38", "mm10", or "mm39"
#' @param project_dir Character: Path to project directory (default: current directory)
#' @param force Logical: Overwrite existing .locusPackRat directory if it exists
#'
#' @return Invisible TRUE on success
#'
#' @note Genome coordinates and orthology mappings are loaded from cached data files
#'       in tests/test_data/ (or inst/extdata/ when installed as a package).
#'       Currently supports mm39 (mouse) and hg38 (human) genomes.
#'
#' @importFrom utils head packageVersion
#' @importFrom data.table setnames fread fwrite setDT := copy
#' @importFrom dplyr bind_rows
#' @importFrom jsonlite write_json read_json
#'
#' @export
initPackRat <- function(data,
                       mode = c("gene", "region"),
                       species = c("human", "mouse"),
                       genome = c("hg38", "hg19", "mm39", "mm10"),
                       project_dir = ".",
                       force = FALSE) {

  # Validate inputs
  mode <- match.arg(mode)
  species <- match.arg(species)
  genome <- match.arg(genome)

  # Check species/genome compatibility
  if ((species == "human" && genome %in% c("mm10", "mm39")) ||
      (species == "mouse" && genome %in% c("hg19", "hg38"))) {
    stop("Incompatible species and genome specified")
  }

  # Set up directory structure
  packrat_dir <- file.path(project_dir, ".locusPackRat")

  if (dir.exists(packrat_dir) && !force) {
    stop(".locusPackRat directory already exists. Use force=TRUE to overwrite.")
  }

  message("Initializing locusPackRat project...")

  # Create directory structure
  dirs <- c(
    packrat_dir,
    file.path(packrat_dir, "input"),
    file.path(packrat_dir, "supplementary"),
    file.path(packrat_dir, "cache"),
    file.path(packrat_dir, "output")
  )

  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
    }
  }

  # Convert input to data.table
  setDT(data)

  # Process based on mode
  if (mode == "gene") {
    message("Processing gene list...")
    processed_data <- .processGeneInput(data, species, genome)
    output_file <- file.path(packrat_dir, "input", "genes.csv")
  } else {
    message("Processing region list...")
    processed_data <- .processRegionInput(data, species, genome)
    output_file <- file.path(packrat_dir, "input", "regions.csv")
  }

  # Generate orthology table
  message("Generating orthology information...")
  orthology_data <- .generateOrthologyTable(processed_data, species, genome)

  # Save processed data
  fwrite(processed_data, output_file)
  message(sprintf("Saved %s data to %s", mode, output_file))

  # Save orthology data
  if (!is.null(orthology_data) && nrow(orthology_data) > 0) {
    ortho_file <- file.path(packrat_dir, "input", "orthology.csv")
    fwrite(orthology_data, ortho_file)
    message(sprintf("Saved orthology data to %s", ortho_file))
  }

  # Create config file
  # Get package version safely (may not be installed as package yet)
  pkg_version <- tryCatch(
    packageVersion("locusPackRat"),
    error = function(e) "0.1.0-dev"
  )

  config <- list(
    mode = mode,
    species = species,
    genome = genome,
    initialization_date = Sys.Date(),
    package_version = as.character(pkg_version),
    input_columns = names(processed_data),
    n_entries = nrow(processed_data)
  )

  config_file <- file.path(packrat_dir, "config.json")
  jsonlite::write_json(config, config_file, pretty = TRUE, auto_unbox = TRUE)
  message(sprintf("Created config file: %s", config_file))

  message("\nlocusPackRat project initialized successfully!")
  message(sprintf("Mode: %s | Species: %s | Genome: %s", mode, species, genome))
  message(sprintf("Processed %d %s", nrow(processed_data), ifelse(mode == "gene", "genes", "regions")))

  invisible(TRUE)
} 

#' Add Supplementary Table to locusPackRat Project
#'
#' Adds supplementary data tables linked to genes or regions in an existing project.
#' Can handle gene-based data, genomic regions, or point coordinates (e.g., SNPs).
#'
#' @param data data.frame or data.table with supplementary information
#' @param table_name Character: Name for the output table (without .csv extension)
#' @param link_type Character: "gene", "region", or "point".
#'                  "point" treats data as genomic positions (SNPs) and maps them to overlapping regions.
#' @param link_by Character: Column name in data to use for linking.
#' @param project_dir Character: Path to project directory containing .locusPackRat
#'
#' @return Invisible TRUE on success
#'
#' @importFrom data.table fread fwrite setDT merge.data.table := setnames
#' @importFrom jsonlite read_json write_json
#'
#' @export
addRatTable <- function(data,
                        table_name,
                        link_type = c("gene", "region", "point"),
                        link_by = NULL,
                        project_dir = ".") {

  # Validate inputs
  link_type <- match.arg(link_type)
  # Preserve the original type for the config file
  original_link_type <- link_type
  
  if (missing(table_name) || is.null(table_name)) {
    stop("table_name must be provided")
  }

  # Check for existing project
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    stop("No .locusPackRat directory found. Run initPackRat() first.")
  }

  # Read config
  config_file <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_file)) {
    stop("Config file not found. Project may be corrupted.")
  }
  config <- jsonlite::read_json(config_file)

  message(sprintf("Adding supplementary table to %s %s project...",
                  config$species, config$genome))

  # Convert input to data.table
  setDT(data)

  # --- HANDLE POINT DATA ---
  # If "point", normalize to "region" format (start == end)
  if (link_type == "point") {
    message("Normalizing point data to genomic intervals...")
    
    # 1. Standardize Chromosome column
    if (!"chr" %in% names(data)) {
       # Try standard alternatives (e.g., chrom, chromosome, seqname)
       chr_col <- grep("^chr|^chrom|^seqname", names(data), value=TRUE, ignore.case=TRUE)[1]
       if (!is.na(chr_col)) {
         setnames(data, chr_col, "chr")
         message(sprintf("  Mapping '%s' to 'chr'", chr_col))
       } else {
         stop("Link type 'point' requires a chromosome column (e.g., 'chr', 'chrom')")
       }
    }
    
    # 2. Standardize Position column
    if (!"pos" %in% names(data)) {
       # Try standard alternatives (e.g., pos, bp, location, start)
       pos_col <- grep("^pos|^bp|^loc|^start", names(data), value=TRUE, ignore.case=TRUE)[1]
       if (!is.na(pos_col)) {
         setnames(data, pos_col, "pos")
         message(sprintf("  Mapping '%s' to 'pos'", pos_col))
       } else {
         stop("Link type 'point' requires a position column (e.g., 'pos', 'bp')")
       }
    }
    
    # 3. Create start/end for overlap joining (Point is an interval of width 0 or 1)
    data[, start := as.numeric(pos)]
    data[, end := as.numeric(pos)]
    
    # Switch processing mode to "region" for the rest of the function
    link_type <- "region"
  }

  # Load existing data based on link type and config mode
  if (link_type == "gene") {
    # Check if we have genes to link to
    genes_file <- file.path(packrat_dir, "input", "genes.csv")
    if (config$mode == "region") {
      # For region mode, we need to extract genes from regions first
      regions_file <- file.path(packrat_dir, "input", "regions.csv")
      if (!file.exists(regions_file)) {
        stop("Regions file not found, cannot link gene data")
      }
      existing_data <- fread(regions_file)
      # Get unique genes from regions
      if ("genes" %in% names(existing_data)) {
        # Assuming genes are stored as comma-separated in a column
        all_genes <- unlist(strsplit(existing_data$genes, ","))
        all_genes <- unique(trimws(all_genes))
        existing_data <- data.table(gene_symbol = all_genes)
      } else {
        stop("No gene information found in regions data")
      }
    } else {
      # Gene mode - direct load
      if (!file.exists(genes_file)) {
        stop("Genes file not found")
      }
      existing_data <- fread(genes_file)
    }

    # Perform linking
    message(sprintf("Linking data by %s...", ifelse(is.null(link_by), "auto-detection", link_by)))
    linked_data <- .linkGeneData(data, existing_data, link_by)

  } else {
    linked_data <- data
  }

  # Save supplementary table
  output_file <- file.path(packrat_dir, "supplementary", paste0(table_name, ".csv"))
  fwrite(linked_data, output_file)
  message(sprintf("Saved supplementary table to %s", output_file))
  message(sprintf("Linked %d of %d input rows", nrow(linked_data), nrow(data)))

  # Update config with table info
  if (is.null(config$supplementary_tables)) {
    config$supplementary_tables <- list()
  }

  config$supplementary_tables[[table_name]] <- list(
    link_type = original_link_type, # Store what the user actually requested
    link_by = ifelse(is.null(link_by), "auto", link_by),
    date_added = Sys.Date(),
    n_rows = nrow(linked_data),
    n_cols = ncol(linked_data),
    columns = names(linked_data)
  )

  jsonlite::write_json(config, config_file, pretty = TRUE, auto_unbox = TRUE)
  message("Updated config file")

  invisible(TRUE)
}


#' List Supplementary Tables in locusPackRat Project
#'
#' Returns information about available supplementary tables in a project
#'
#' @param project_dir Character: Path to project directory containing .locusPackRat
#'
#' @return data.table with table information (name, link_type, link_by, n_rows, n_cols)
#'
#' @importFrom data.table data.table fread
#' @importFrom dplyr bind_rows
#' @importFrom jsonlite read_json
#'
#' @export
listPackRatTables <- function(project_dir = ".") {
  # Check for existing project
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    stop("No .locusPackRat directory found. Run initPackRat() first.")
  }

  # Read config
  config_file <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_file)) {
    warning("Config file not found. Project may be corrupted.")
    config <- list()
  } else {
    config <- jsonlite::read_json(config_file)
  }

  # Check supplementary directory
  supp_dir <- file.path(packrat_dir, "supplementary")
  if (!dir.exists(supp_dir)) {
    message("No supplementary tables found")
    return(data.table())
  }

  # Get list of tables
  supp_files <- list.files(supp_dir, pattern = "\\.csv$", full.names = FALSE)
  if (length(supp_files) == 0) {
    message("No supplementary tables found")
    return(data.table())
  }

  # Build table information
  table_info <- data.table()
  for (sf in supp_files) {
    table_name <- sub("\\.csv$", "", sf)

    # Get info from config if available
    if (!is.null(config$supplementary_tables) &&
        table_name %in% names(config$supplementary_tables)) {
      tbl_config <- config$supplementary_tables[[table_name]]
      info <- data.table(
        table_name = table_name,
        link_type = ifelse(is.null(tbl_config$link_type), NA_character_, tbl_config$link_type),
        link_by = ifelse(is.null(tbl_config$link_by), NA_character_, tbl_config$link_by),
        n_rows = ifelse(is.null(tbl_config$n_rows), NA_integer_, tbl_config$n_rows),
        n_cols = ifelse(is.null(tbl_config$n_cols), NA_integer_, tbl_config$n_cols),
        date_added = as.character(ifelse(is.null(tbl_config$date_added), NA_character_, tbl_config$date_added))
      )
    } else {
      # Try to get basic info from file
      file_path <- file.path(supp_dir, sf)
      tryCatch({
        temp_data <- fread(file_path, nrows = 1)
        info <- data.table(
          table_name = table_name,
          link_type = NA_character_,
          link_by = NA_character_,
          n_rows = as.integer(system(sprintf("wc -l < '%s'", file_path), intern = TRUE)) - 1L,
          n_cols = ncol(temp_data),
          date_added = NA_character_
        )
      }, error = function(e) {
        info <- data.table(
          table_name = table_name,
          link_type = NA_character_,
          link_by = NA_character_,
          n_rows = NA_integer_,
          n_cols = NA_integer_,
          date_added = NA_character_
        )
      })
    }

    table_info <- bind_rows(table_info, info)
  }

  message(sprintf("Found %d supplementary table(s):", nrow(table_info)))
  for (i in seq_len(nrow(table_info))) {
    message(sprintf("  - %s: %d rows with %d cols, linked by '%s'",
                   table_info$table_name[i],
                   table_info$n_rows[i],
                   table_info$n_cols[i],
                   table_info$link_by[i]))
  }

  return(table_info)
}

#' Collapse Column Values for Merging
#' @noRd
.collapse_column <- function(x) {
  # Drop NAs and empty strings (if character)
  if (is.character(x)) {
    x <- x[x != ""]
  }
  x <- unique(x[!is.na(x)])

  # If no values left, return NA of the *same type* as the column
  if (length(x) == 0L) {
    return(x[NA_integer_])   # NA with same typeof(x)
  }

  # If only one value, just return it (keeps original type)
  if (length(x) == 1L) {
    return(x)
  }

  # If character and multiple values, paste them
  if (is.character(x)) {
    return(paste(x, collapse = "; "))
  }

  # For non-character with multiple distinct values, just keep the first
  # (keeps type numeric/logical/etc, avoids mixing with character)
  x[1L]
}

#' Generate Gene Sheet from locusPackRat Project
#'
#' Creates formatted gene sheets (CSV or Excel) from project data with optional filtering and highlighting
#'
#' @param filter_expr Character: Expression to filter genes (uses data.table syntax, e.g., "chr == 1 & start > 1000000")
#' @param highlight_genes Character vector: Gene symbols to highlight
#' @param highlight_expr Character: Expression to determine which genes to highlight
#' @param highlight_color Character: Background color for highlighted genes (hex color)
#' @param summary_level Character: "gene" for one row per gene, "locus" for one row per locus/region
#' @param output_file Character: Output file path (defaults to .locusPackRat/output/gene_sheet.csv or .xlsx)
#' @param format Character: "csv" for plain CSV, "excel" for formatted Excel workbook
#' @param split_by Character: For Excel, how to split into tabs - "none", "chromosome", "criteria"
#' @param split_criteria Named list: For split_by="criteria", list of tab_name = filter_expression
#' @param include_supplementary Logical or character vector: TRUE to include all supplementary tables,
#'                               FALSE for none, or character vector of specific table names to include
#' @param collapse_multiples Logical: when merging supplementary tables,
#'   collapse multiple matches per key into a single row by aggregating
#'   values (e.g. concatenating text). If FALSE, keep all matches as
#'   separate rows.
#' @param project_dir Character: Path to project directory containing .locusPackRat
#'
#' @return data.table of the gene sheet (invisible)
#'
#' @importFrom data.table fread fwrite setDT merge.data.table := .SD
#' @importFrom jsonlite read_json
#' @importFrom openxlsx createWorkbook addWorksheet writeDataTable createStyle addStyle setColWidths freezePane saveWorkbook
#'
#' @export
makeGeneSheet <- function(filter_expr = NULL,
                          highlight_genes = NULL,
                          highlight_expr = NULL,
                          highlight_color = "#FFFF99",
                          summary_level = c("gene", "locus"),
                          output_file = NULL,
                          format = c("excel", "csv"),
                          split_by = c("none", "chromosome", "criteria"),
                          split_criteria = NULL,
                          include_supplementary = TRUE,
                          collapse_multiples = TRUE,
                          project_dir = ".") {

  summary_level <- match.arg(summary_level)
  format <- match.arg(format)
  split_by <- match.arg(split_by)

  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) stop(sprintf("No .locusPackRat directory found at %s", packrat_dir))

  config <- jsonlite::read_json(file.path(packrat_dir, "config.json"))
  message(sprintf("Generating %s sheet from %s...", config$mode, packrat_dir))

  # 1. LOAD BASE DATA
  if (config$mode == "gene") {
    base_data <- fread(file.path(packrat_dir, "input", "genes.csv"))
    if(!"pk_id" %in% names(base_data)) base_data[, pk_id := .I]
  } else {
    region_data <- fread(file.path(packrat_dir, "input", "regions.csv"))
    base_data <- .expandRegionsToGenes(region_data)
    if(!"pk_id" %in% names(base_data)) base_data[, pk_id := .I]
  }

  # 2. LOAD ORTHOLOGY
  ortho_file <- file.path(packrat_dir, "input", "orthology.csv")
  if (file.exists(ortho_file)) {
    ortho_data <- fread(ortho_file)
    if ("ortho_species" %in% names(ortho_data)) ortho_data[, ortho_species := NULL]

    if (config$species == "mouse") {
      mx <- "gene_symbol"; my <- "mouse_gene_symbol"; keep <- c("human_gene_symbol", "human_ensembl_id")
    } else {
      mx <- "gene_symbol"; my <- "human_gene_symbol"; keep <- c("mouse_gene_symbol", "mouse_ensembl_id")
    }

    if (mx %in% names(base_data) && my %in% names(ortho_data)) {
      base_data <- merge(base_data, ortho_data[, c(my, keep), with=FALSE], by.x = mx, by.y = my, all.x = TRUE)
    }
  }

  # 3. LINK SUPPLEMENTARY TABLES
  if (!isFALSE(include_supplementary)) {
    supp_dir <- file.path(packrat_dir, "supplementary")
    
    if (dir.exists(supp_dir)) {
      all_supp <- list.files(supp_dir, pattern = "\\.csv$", full.names = FALSE)
      
      # Debug: Check what files are found
      if(length(all_supp) == 0) {
        warning(sprintf("No CSV files found in supplementary directory: %s", supp_dir))
      } else {
        message(sprintf("Found %d supplementary files: %s", length(all_supp), paste(all_supp, collapse=", ")))
      }

      if (is.character(include_supplementary)) {
        tables_to_include <- intersect(paste0(include_supplementary, ".csv"), all_supp)
      } else {
        tables_to_include <- all_supp
      }

      for (sf in tables_to_include) {
        table_name <- sub("\\.csv$", "", sf)
        supp_data <- fread(file.path(supp_dir, sf))
        
        # --- A. DETECT KEYS ---
        merge_cols <- NULL
        # Check Config
        if (!is.null(config$supplementary_tables[[table_name]]$link_by) && config$supplementary_tables[[table_name]]$link_by != "auto") {
          merge_cols <- trimws(strsplit(config$supplementary_tables[[table_name]]$link_by, ",")[[1]])
        }
        # Fallback Auto-detect
        if (is.null(merge_cols)) {
          common <- intersect(names(base_data), names(supp_data))
          if ("region_id" %in% common) merge_cols <- "region_id"
          else if (all(c("chr", "start", "end") %in% common)) merge_cols <- c("chr", "start", "end")
          else if ("gene_symbol" %in% common) merge_cols <- "gene_symbol"
        }

        # --- B. EXECUTE MERGE LOGIC ---
        if (!is.null(merge_cols)) {
          # Capture columns we EXPECT to add (for safety filling later)
          expected_cols <- setdiff(names(supp_data), merge_cols)

          # CASE 1: COORDINATE OVERLAP (Interval Join)
          if (all(c("chr", "start", "end") %in% merge_cols)) {
             message(sprintf("  Linking %s using Interval Overlap...", table_name))
             
             # Normalize Coordinates
             supp_data[, chr := gsub("^chr", "", as.character(chr), ignore.case=TRUE)]
             supp_data[, start := as.integer(start)][, end := as.integer(end)]
             base_data[, chr := gsub("^chr", "", as.character(chr), ignore.case=TRUE)]
             base_data[, start := as.integer(start)][, end := as.integer(end)]

             # Non-Equi Join
             # Note: This join might produce 0 rows if no overlaps
             overlap_join <- supp_data[base_data, 
                                       on = .(chr, start <= end, end >= start), 
                                       nomatch = NULL] 
             
             val_cols <- setdiff(names(supp_data), c("chr", "start", "end"))
             
             if (length(val_cols) > 0) {
               # Aggregation
              collapsed_vals <- overlap_join[
                ,
                lapply(.SD, .collapse_column),
                by = pk_id,
                .SDcols = val_cols
              ]
               
               # Merge
               base_data <- merge(base_data, collapsed_vals, by = "pk_id", all.x = TRUE)
             }

          } else {
            # CASE 2: EXACT MERGE
            message(sprintf("  Linking %s (Exact Match: %s)", table_name, paste(merge_cols, collapse=",")))
            
            if (collapse_multiples) {
              val_cols <- setdiff(names(supp_data), merge_cols)
              if (length(val_cols) > 0) {
                supp_data <- supp_data[
                  ,
                  lapply(.SD, .collapse_column),
                  by = merge_cols,
                  .SDcols = val_cols
                ]
              } else {
                supp_data <- unique(supp_data)
              }
            }
            base_data <- merge(base_data, supp_data, by = merge_cols, all.x = TRUE, suffixes = c("", paste0(".", table_name)))
          }
          
          # SAFETY CHECK: Ensure expected columns exist
          # If merge failed or yielded 0 matches, columns might be missing.
          # We force add them as NA to prevent 'Object not found' errors in filtering.
          missing_cols <- setdiff(expected_cols, names(base_data))
          if(length(missing_cols) > 0) {
             message(sprintf("    Warning: %s matched 0 rows. Adding empty columns: %s", table_name, paste(missing_cols, collapse=", ")))
             base_data[, (missing_cols) := NA]
          }

          # Clean suffixes
          bad_cols <- grep(paste0("\\.", table_name, "$"), names(base_data), value = TRUE)
          if(length(bad_cols) > 0) base_data[, (bad_cols) := NULL]

        } else {
          warning(sprintf("  Skipping %s: No link columns found.", table_name))
        }
      } 
    }
  }

  if("pk_id" %in% names(base_data)) base_data[, pk_id := NULL]

  # 4. FILTERING
  if (!is.null(filter_expr)) {
    tryCatch({
      base_data <- base_data[eval(parse(text = filter_expr))]
    }, error = function(e) stop("Filter Expression Error: ", e$message))
  }

  # 5. SUMMARIZE & SAVE
  output_data <- if (summary_level == "locus") .summarizeAtLocusLevel(base_data, config) else base_data

  output_dir <- file.path(packrat_dir, "output")
  if (is.null(output_file)) output_file <- paste0("gene_sheet", ifelse(format == "excel", ".xlsx", ".csv"))
  final_path <- file.path(output_dir, basename(output_file))

  if (format == "csv") {
    fwrite(output_data, final_path)
    message(sprintf("Saved CSV to %s", final_path))
  } else {
    wb <- .createFormattedExcel(output_data, highlight_genes, highlight_expr, highlight_color, split_by, split_criteria)
    openxlsx::saveWorkbook(wb, final_path, overwrite = TRUE)
    message(sprintf("Saved Excel to %s", final_path))
  }

  invisible(output_data)
}

# ========== Helper Functions ==========

#' Process Gene Input Data
#' @noRd
.processGeneInput <- function(data, species, genome) {
  # Detect ID type based on column names
  gene_cols <- c("gene_symbol", "gene_name", "symbol", "ensembl_id", "gene_id", "id")
  id_col <- intersect(names(data), gene_cols)[1]

  if (is.na(id_col)) {
    stop("Could not find gene identifier column. Expected one of: ",
         paste(gene_cols, collapse = ", "))
  }

  # Standardize column names
  setnames(data, id_col, "input_id", skip_absent = TRUE)

  # Load genome coordinate data from cached files
  coords_file <- NULL
  if (species == "mouse") {
    if (genome == "mm39") {
      coords_file <- system.file("extdata", "mouse_coords_mm39.csv", package = "locusPackRat")
      if (coords_file == "") {
        # Try test data location if package not installed
        coords_file <- "tests/test_data/mouse_coords_mm39.csv"
      }
    }
    # Could add mm10 support here
  } else if (species == "human") {
    if (genome == "hg38") {
      coords_file <- system.file("extdata", "human_coords_hg38.csv", package = "locusPackRat")
      if (coords_file == "") {
        # Try test data location if package not installed
        coords_file <- "tests/test_data/human_coords_hg38.csv"
      }
    }
    # Could add hg19 support here
  }

  # Initialize result with input gene symbols
  result <- data.table(
    gene_symbol = data$input_id,
    ensembl_id = NA_character_,
    chr = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    strand = NA_character_
  )

  # Load coordinates if file exists
  if (!is.null(coords_file) && file.exists(coords_file)) {
    coords_data <- fread(coords_file)

    # Match on gene_symbol or ensembl_id
    id_type <- ifelse(grepl("^ENS", data$input_id[1]), "ensembl_id", "gene_symbol")

    if (id_type == "gene_symbol") {
      # Merge by gene symbol
      result <- merge(result[, .(gene_symbol)], coords_data,
                     by = "gene_symbol", all.x = TRUE)
    } else {
      # Merge by ensembl ID
      result[, ensembl_id := data$input_id]
      result <- merge(result[, .(ensembl_id)], coords_data,
                     by = "ensembl_id", all.x = TRUE)
    }

    # Ensure numeric coordinates
    if ("start" %in% names(result)) result[, start := as.integer(start)]
    if ("end" %in% names(result)) result[, end := as.integer(end)]
  } else {
    message("  Note: Coordinate data file not found, using NA placeholders")
  }

  # Preserve any additional columns from input data (but not the original id column)
  extra_cols <- setdiff(names(data), c(id_col, "input_id"))
  if (length(extra_cols) > 0) {
    # Merge back the extra columns
    data_with_extras <- copy(data)
    data_with_extras[, gene_symbol := input_id]
    result <- merge(result, data_with_extras[, c("gene_symbol", extra_cols), with = FALSE],
                   by = "gene_symbol", all.x = TRUE)
  }

  return(result)
}

#' Process Region Input Data
#' @noRd
.processRegionInput <- function(data, species, genome) {
  # Check for required columns
  required <- c("chr", "start", "end")
  missing <- setdiff(required, names(data))

  if (length(missing) > 0) {
    # Try alternative names
    alt_names <- list(
      chr = c("chrom", "chromosome", "seqname"),
      start = c("start_pos", "chromStart", "begin"),
      end = c("end_pos", "chromEnd", "stop")
    )

    for (req in missing) {
      alts <- alt_names[[req]]
      found <- intersect(names(data), alts)
      if (length(found) > 0) {
        setnames(data, found[1], req)
      }
    }
  }
  # Final check
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns for regions: ", paste(missing, collapse = ", "))
  }

  # Ensure numeric coordinates and character chromosome
  data[, start := as.numeric(start)]
  data[, end := as.numeric(end)]
  data[, chr := as.character(chr)]

  # Add region ID if not present
  if (!"region_id" %in% names(data)) {
    data[, region_id := paste0("region_", .I)]
  }
  # Get genes in regions
  coord_type <- paste0(species, "_coords_", genome, ".csv")
  coords_file <- system.file("extdata", coord_type, package = "locusPackRat")
  coords <- fread(coords_file)

  # Ensure numeric coordinates
  coords[, start := as.numeric(start)]
  coords[, end := as.numeric(end)]
  coords[, chr := as.character(chr)]
  
  # look for overlaps
  overlaps <- coords[data,
                      .(region_id = i.region_id, gene_symbol = x.gene_symbol), 
                      on = .(chr, start <= end, end >= start),
                      nomatch = NULL]
  if (nrow(overlaps) > 0) {
      # Collapse multiple genes into a single comma-separated string
      genes_aggregated <- overlaps[, .(genes = paste(unique(gene_symbol), collapse = ", ")),
                                  by = region_id]
      
      # Merge back to the original data
      data <- merge(data, genes_aggregated, by = "region_id", all.x = TRUE)
      
  } else {
      # If no overlaps found, create empty column
      data[, genes := NA_character_]
  }
return(data)
}

#' Generate Orthology Table
#' @noRd
.generateOrthologyTable <- function(data, species, genome) {
  # Load orthology data from cached file
  ortho_file <- system.file("extdata", "orthology.csv", package = "locusPackRat")
  if (ortho_file == "") {
    # Try test data location if package not installed
    ortho_file <- "tests/test_data/orthology_mapping.csv"
  }

  if (!file.exists(ortho_file)) {
    message("  Note: Orthology data file not found, using NA placeholders")

    # Fallback to placeholder data
    if ("gene_symbol" %in% names(data)) {
      genes <- unique(data$gene_symbol)
      genes <- genes[!is.na(genes)]

      if (length(genes) > 0) {
        ortho <- data.table(
          source_gene_symbol = genes,
          source_ensembl_id = NA_character_,
          ortho_gene_symbol = NA_character_,
          ortho_ensembl_id = NA_character_,
          ortho_species = ifelse(species == "human", "mouse", "human")
        )

        # Add species-specific column names
        if (species == "human") {
          setnames(ortho,
                  c("source_gene_symbol", "source_ensembl_id", "ortho_gene_symbol", "ortho_ensembl_id"),
                  c("human_gene_symbol", "human_ensembl_id", "mouse_gene_symbol", "mouse_ensembl_id"))
        } else {
          setnames(ortho,
                  c("source_gene_symbol", "source_ensembl_id", "ortho_gene_symbol", "ortho_ensembl_id"),
                  c("mouse_gene_symbol", "mouse_ensembl_id", "human_gene_symbol", "human_ensembl_id"))
        }

        return(ortho)
      }
    }
    return(NULL)
  }

  # Load orthology mapping data
  ortho_mapping <- fread(ortho_file)

  # Remove duplicates if any
  ortho_mapping <- unique(ortho_mapping)

  # Get unique genes from input data
  if ("gene_symbol" %in% names(data)) {
    genes <- unique(data$gene_symbol)
    genes <- genes[!is.na(genes)]
  } else if ("ensembl_id" %in% names(data)) {
    ensembl_ids <- unique(data$ensembl_id)
    ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]
    # We'll match on Ensembl IDs instead
    genes <- NULL
  } else {
    return(NULL)
  }

  # Filter orthology data for our genes
  if (species == "mouse") {
    # Filter for mouse genes
    if (!is.null(genes)) {
      result <- ortho_mapping[mouse_gene_symbol %in% genes]
    } else {
      result <- ortho_mapping[mouse_ensembl_id %in% ensembl_ids]
    }

    # Ensure we have all input genes in result (even if no ortholog)
    if (!is.null(genes)) {
      missing_genes <- setdiff(genes, result$mouse_gene_symbol)
      if (length(missing_genes) > 0) {
        missing_rows <- data.table(
          mouse_gene_symbol = missing_genes,
          mouse_ensembl_id = NA_character_,
          human_gene_symbol = NA_character_,
          human_ensembl_id = NA_character_
        )
        result <- bind_rows(result, missing_rows)
      }
    }
  } else {
    # Filter for human genes
    if (!is.null(genes)) {
      result <- ortho_mapping[human_gene_symbol %in% genes]
    } else {
      result <- ortho_mapping[human_ensembl_id %in% ensembl_ids]
    }

    # Ensure we have all input genes in result (even if no ortholog)
    if (!is.null(genes)) {
      missing_genes <- setdiff(genes, result$human_gene_symbol)
      if (length(missing_genes) > 0) {
        missing_rows <- data.table(
          human_gene_symbol = missing_genes,
          human_ensembl_id = NA_character_,
          mouse_gene_symbol = NA_character_,
          mouse_ensembl_id = NA_character_
        )
        result <- bind_rows(result, missing_rows)
      }
    }
  }

  # Remove the ortho_species column if it exists (it's redundant)
  if ("ortho_species" %in% names(result)) {
    result[, ortho_species := NULL]
  }

  return(result)
}

#' Link Gene Data
#' @noRd
.linkGeneData <- function(new_data, existing_genes, link_by = NULL) {
  if (is.null(link_by)) {
    # Auto-detect linking column
    possible_links <- c("gene_symbol", "ensembl_id", "gene_name", "gene_id")
    link_cols <- intersect(names(new_data), possible_links)

    if (length(link_cols) == 0) {
      stop("Could not auto-detect linking column. Please specify link_by parameter.")
    }
    link_by <- link_cols[1]
  }

  # Check if link column exists in both datasets
  if (!link_by %in% names(new_data)) {
    stop(sprintf("Link column '%s' not found in input data", link_by))
  }

  # Find matching column in existing data
  existing_link <- intersect(names(existing_genes), c(link_by, "gene_symbol", "ensembl_id"))
  if (length(existing_link) == 0) {
    stop("Could not find matching column in existing gene data")
  }

  # Perform merge to find matching rows
  # But only keep the linking column and new data columns
  matched <- merge(existing_genes[, ..existing_link], new_data,
                 by.x = existing_link[1], by.y = link_by,
                 all = FALSE)

  # Rename the linking column back to its original name if needed
  if (existing_link[1] != link_by) {
    setnames(matched, existing_link[1], link_by)
  }

  return(matched)
}

#' Link Region Data
#' @importFrom data.table foverlaps setkey
#' @noRd
.linkRegionData <- function(new_data, existing_regions, genome) {
  # Check for coordinate columns in new data
  coord_cols <- c("chr", "start", "end")
  if (!all(coord_cols %in% names(new_data))) {
    stop("Region linking requires chr, start, end columns in input data")
  }

  # Convert to numeric
  new_data[, start := as.numeric(start)]
  new_data[, end := as.numeric(end)]
  existing_regions[, start := as.numeric(start)]
  existing_regions[, end := as.numeric(end)]

  # Force chromosome to character to avoid integer/character join errors
  # (fread often guesses integer for simple chromosomes like "1", "2")
  new_data[, chr := as.character(chr)]
  existing_regions[, chr := as.character(chr)]

  # Find overlaps
  setkey(new_data, chr, start, end)
  setkey(existing_regions, chr, start, end)

  # Use data.table's foverlaps for efficient overlap detection
  overlaps <- data.table::foverlaps(new_data, existing_regions, nomatch = NULL)

  return(overlaps)
}

#' Expand Regions to Genes
#' @noRd
.expandRegionsToGenes <- function(region_data) {
  if (!"genes" %in% names(region_data) || all(is.na(region_data$genes)) || all(region_data$genes == "")) {
    warning("`genes` column is missing or empty, returning empty gene table")
    return(data.table())
  }

  region_data[, genes := as.character(genes)]
  # Split comma-separated genes
  genes_list <- strsplit(region_data$genes, ", ")

  # Create expanded table
  result <- region_data[rep(seq_len(nrow(region_data)), lengths(genes_list))]
  result[, gene_symbol := unlist(genes_list)]
  result[, gene_symbol := trimws(gene_symbol)]

  # Add region info
  result[, from_region := paste0(chr, ":", start, "-", end)]

  return(result)
}

#' Summarize at Locus Level
#' @noRd
.summarizeAtLocusLevel <- function(data, config) {
  if (config$mode == "region") {
    # Group by region
    if ("region_id" %in% names(data)) {
      summary <- data[, .(
        n_genes = .N,
        genes = paste(unique(gene_symbol), collapse = ", "),
        chr = chr[1],
        start = min(start, na.rm = TRUE),
        end = max(end, na.rm = TRUE)
      ), by = region_id]
    } else {
      # Group by coordinates
      summary <- data[, .(
        n_genes = .N,
        genes = paste(unique(gene_symbol), collapse = ", ")
      ), by = .(chr, start, end)]
    }
  } else {
    # For gene mode, group by chromosome
    summary <- data[, .(
      n_genes = .N,
      genes = paste(unique(gene_symbol), collapse = ", "),
      start = min(start, na.rm = TRUE),
      end = max(end, na.rm = TRUE)
    ), by = chr]
  }

  return(summary)
}

#' Create Formatted Excel Workbook (Safely)
#' @noRd
.createFormattedExcel <- function(data, highlight_genes, highlight_expr,
                                  highlight_color, split_by, split_criteria) {
  
  wb <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(textDecoration = "bold", halign = "center", border = "bottom", borderStyle = "medium")
  highlight_style <- openxlsx::createStyle(fgFill = highlight_color)
  stripe_style <- openxlsx::createStyle(fgFill = "#F0F0F0")
  
  # Determine sheets
  sheets <- list()
  
  if (split_by == "none") {
    sheets[["All Genes"]] <- data
  } else if (split_by == "chromosome") {
    if ("chr" %in% names(data)) {
      sheets <- split(data, data$chr)
      names(sheets) <- paste0("Chr_", names(sheets))
    } else {
      sheets[["All Genes"]] <- data
    }
  } else if (split_by == "criteria" && !is.null(split_criteria)) {
    sheets[["All_Data"]] <- data
    for (tab_name in names(split_criteria)) {
      f_expr <- split_criteria[[tab_name]]
      tryCatch({
        # Verify columns exist before eval
        # Note: eval might still fail if logic is wrong, but Object Not Found is usually missing cols.
        filtered <- data[eval(parse(text = f_expr))]
        sheets[[tab_name]] <- filtered
      }, error = function(e) {
        # Detailed warning
        warning(sprintf("Skipping tab '%s': %s", tab_name, e$message))
      })
    }
  }
  
  for (sheet_name in names(sheets)) {
    sheet_data <- sheets[[sheet_name]]
    if (nrow(sheet_data) == 0) next
    
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeDataTable(wb, sheet = sheet_name, x = sheet_data, tableName = gsub("[^A-Za-z0-9]", "_", sheet_name), withFilter = TRUE)
    openxlsx::addStyle(wb, sheet = sheet_name, style = header_style, rows = 1, cols = 1:ncol(sheet_data))
    
    if (nrow(sheet_data) > 0) {
      rows <- seq(2, nrow(sheet_data) + 1, by = 2)
      openxlsx::addStyle(wb, sheet = sheet_name, style = stripe_style, rows = rows, cols = 1:ncol(sheet_data), gridExpand = TRUE)
    }
    
    rows_to_hl <- c()
    if (!is.null(highlight_genes) && "gene_symbol" %in% names(sheet_data)) {
      rows_to_hl <- which(sheet_data$gene_symbol %in% highlight_genes) + 1
    }
    if (!is.null(highlight_expr)) {
      try({
        expr_rows <- which(sheet_data[, eval(parse(text = highlight_expr))]) + 1
        rows_to_hl <- unique(c(rows_to_hl, expr_rows))
      }, silent = TRUE)
    }
    if (length(rows_to_hl) > 0) {
      openxlsx::addStyle(wb, sheet = sheet_name, style = highlight_style, rows = rows_to_hl, cols = 1:ncol(sheet_data), gridExpand = TRUE)
    }
    
    openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:ncol(sheet_data), widths = "auto")
    openxlsx::freezePane(wb, sheet = sheet_name, firstActiveRow = 2)
  }
  return(wb)
}