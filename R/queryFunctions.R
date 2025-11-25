#' @noRd
.getProjectGenes <- function(packrat_dir) {
  config <- jsonlite::read_json(file.path(packrat_dir, "config.json"))
  
  if (config$mode == "gene") {
    dt <- data.table::fread(file.path(packrat_dir, "input", "genes.csv"))
    return(unique(dt$gene_symbol))
  } else {
    # In region mode, genes are nested in the regions file
    dt <- data.table::fread(file.path(packrat_dir, "input", "regions.csv"))
    # Expand comma-separated genes
    genes <- unlist(strsplit(dt$genes, ", "))
    return(unique(trimws(genes)))
  }
}

#' Query MouseMine for Phenotypes
#' 
#' Draws Data on gene annotations from MouseMine.  MouseMine is downloadable using devtools::install_github("intermine/InterMineR")
#'
#' @param project_dir Character path to the project directory
#' @param limit Integer (optional). If provided, limits the query to the first N genes. Useful for testing.
#' @param chunk_size Integer. Number of genes to query per API call. Default 200.
#'
#' @return Invisible TRUE on success
#' 
#' @export
queryMouseMine <- function(project_dir = ".", limit = NULL, chunk_size = 200) {
  
  # 1. Check Dependencies
  if (!requireNamespace("InterMineR", quietly = TRUE)) {
    stop("Please install 'InterMineR' to use this feature.")
  }
  
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) stop("Project not found.")
  
  # 2. Get Genes
  message("Retrieving gene list from project...")
  gene_syms <- .getProjectGenes(packrat_dir)
  
  # Limit for testing
  if (!is.null(limit)) {
    message(sprintf("TEST MODE: Limiting query to first %d genes", limit))
    gene_syms <- head(gene_syms, limit)
  }
  
  # 3. Setup InterMine
  im <- InterMineR::initInterMine(mine = InterMineR::listMines()["MouseMine"])
  queryGenePath <- InterMineR::getTemplateQuery(im = im, name = "_Feature_Phenotype")
  
  # 4. Batch Processing
  chunks <- split(gene_syms, ceiling(seq_along(gene_syms) / chunk_size))
  
  results_list <- list()
  
  message(sprintf("Querying %d genes in %d batches...", length(gene_syms), length(chunks)))
  
  for (i in seq_along(chunks)) {
    syms <- chunks[[i]]
    
    # Update query constraint
    queryGenePath$where[[4]] <- list(
      path = "OntologyAnnotation.subject",
      op = "LOOKUP",
      value = paste(syms, collapse = ","),
      code = "B"
    )
    
    message(sprintf("  Processing batch %d/%d...", i, length(chunks)))
    
    tryCatch({
      res <- InterMineR::runQuery(im, queryGenePath)
      if (nrow(res) > 0) results_list[[i]] <- res
    }, error = function(e) warning("Batch failed: ", e$message))
    
    Sys.sleep(0.2) # Rate limiting
  }
  
  final_dt <- data.table::rbindlist(results_list, fill = TRUE)
  
  # 5. Clean & Save
  if (nrow(final_dt) > 0) {
    data.table::setnames(final_dt, "OntologyAnnotation.subject.symbol", "gene_symbol", skip_absent=TRUE)
    
    addRatTable(
      data = final_dt,
      table_name = "mouse_phenotypes",
      link_type = "gene",
      link_by = "gene_symbol",
      project_dir = project_dir
    )
    message("Success: 'mouse_phenotypes' added to project.")
  } else {
    warning("No phenotypes found for these genes.")
  }
  
  invisible(final_dt)
}

#' Query Open Targets Genetics
#'
#' Fetches disease associations, genetic constraints, and drug tractability
#' information from the Open Targets GraphQL API.
#'
#' @param project_dir Character path to the project directory
#' @param limit Integer (optional). Limit query to first N genes for testing.
#' @param chunk_size Integer. Number of genes per API call (default 50).
#'
#' @return Invisible TRUE on success
#' 
#' @importFrom stats na.omit  
#' @export
queryOpenTargets <- function(project_dir = ".", limit = NULL, chunk_size = 50) {
  
  # 1. Setup and Checks
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for this function.")
  }
  
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) stop("Project not found.")
  
  # 2. Read Config
  config_path <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_path)) stop("config.json not found in project.")
  config <- jsonlite::read_json(config_path)
  
  message(sprintf("Project Config: %s - %s", config$species, config$genome))

  # 3. Get Genes and Resolve Human IDs
  gene_syms <- .getProjectGenes(packrat_dir)
  
  # Prepare mapping table: gene_symbol (Input) -> human_ensembl_id (For API)
  id_map <- data.table::data.table(gene_symbol = gene_syms)
  
  if (config$species == "human") {
    # --- HUMAN ---
    coords_file <- system.file("extdata", "human_coords_hg38.csv", package = "locusPackRat")
    if (coords_file == "") coords_file <- "inst/extdata/human_coords_hg38.csv" # Fallback for dev
    
    if (!file.exists(coords_file)) stop("Internal 'human_coords_hg38.csv' not found.")
    
    ref_coords <- data.table::fread(coords_file)
    id_map <- merge(id_map, ref_coords[, .(gene_symbol, human_ensembl_id = ensembl_id)], 
                    by = "gene_symbol", all.x = TRUE)
    
  } else if (config$species == "mouse") {
    # --- MOUSE ---
    # 1. Load Mouse Coords (Symbol -> ENSMUSG)
    coords_file <- system.file("extdata", "mouse_coords_mm39.csv", package = "locusPackRat")
    if (coords_file == "") coords_file <- "inst/extdata/mouse_coords_mm39.csv"
    if (!file.exists(coords_file)) stop("Internal 'mouse_coords_mm39.csv' not found.")
    
    ref_coords <- data.table::fread(coords_file)
    
    # 2. Load Internal Orthology (ENSMUSG -> ENSG)
    ortho_file <- system.file("extdata", "orthology.csv", package = "locusPackRat")
    if (ortho_file == "") ortho_file <- "inst/extdata/orthology.csv"
    if (!file.exists(ortho_file)) stop("Internal 'orthology.csv' not found.")
    
    ortho_dat <- data.table::fread(ortho_file)
    
    # 3. Double Join: Input -> Mouse ID -> Human ID
    step1 <- merge(id_map, ref_coords[, .(gene_symbol, mouse_ensembl_id = ensembl_id)], 
                   by = "gene_symbol", all.x = TRUE)
    
    final_map <- merge(step1, ortho_dat[, .(mouse_ensembl_id, human_ensembl_id)], 
                       by = "mouse_ensembl_id", all.x = TRUE)
    
    id_map <- final_map[, .(gene_symbol, human_ensembl_id)]
  } else {
    stop("Unsupported species in config: ", config$species)
  }

  query_ids <- unique(na.omit(id_map$human_ensembl_id))
  
  if (length(query_ids) == 0) {
    warning("No valid Human Ensembl IDs found for these genes. Skipping query.")
    return(NULL)
  }

  if (!is.null(limit)) {
    message(sprintf("TEST MODE: Limiting to first %d genes", limit))
    query_ids <- head(query_ids, limit)
  }
  
  # 4. Define GraphQL Query
  build_query <- function(ids) {
    blocks <- vapply(ids, function(id) {
      alias <- paste0("g_", gsub("-", "_", id))
      sprintf(
        '%s: target(ensemblId: "%s") {
          id
          approvedSymbol
          geneticConstraint { constraintType exp obs score oe oeLower oeUpper }
          tractability { label modality value }
          associatedDiseases(page: { size: 50, index: 0 }) { rows { disease { id name } score } }
        }', alias, id)
    }, FUN.VALUE = character(1))
    paste0("query multiTarget {\n", paste(blocks, collapse = "\n"), "\n}")
  }
  
  # 5. Batch Processing
  batches <- split(query_ids, ceiling(seq_along(query_ids) / chunk_size))
  endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  results_assoc <- list()
  results_constr <- list()
  results_tract <- list()
  
  message(sprintf("Querying Open Targets for %d genes in %d batches...", length(query_ids), length(batches)))
  
  for (i in seq_along(batches)) {
    batch_ids <- batches[[i]]
    message(sprintf("  Processing batch %d/%d...", i, length(batches)))
    
    tryCatch({
      response <- httr::POST(endpoint, body = list(query = build_query(batch_ids)), encode = "json")
      httr::stop_for_status(response)
      
      content <- httr::content(response, "text", encoding="UTF-8")
      data_json <- jsonlite::fromJSON(content, flatten=FALSE)$data
      
      for (alias in names(data_json)) {
        g_dat <- data_json[[alias]]
        if (is.null(g_dat)) next
        
        h_id <- g_dat$id
        
        # --- A. Diseases (Fixed Logic) ---
        if (!is.null(g_dat$associatedDiseases$rows) && length(g_dat$associatedDiseases$rows) > 0) {
          # as.data.table automatically flattens 'disease' {id, name} -> 'disease.id', 'disease.name'
          d_dt <- data.table::as.data.table(g_dat$associatedDiseases$rows)
          
          # Rename dot-separated columns if they exist
          if ("disease.id" %in% names(d_dt)) {
            data.table::setnames(d_dt, 
                                 old = c("disease.id", "disease.name"), 
                                 new = c("disease_id", "disease_name"), 
                                 skip_absent = TRUE)
          }
          
          d_dt[, human_ensembl_id := h_id]
          results_assoc[[length(results_assoc) + 1]] <- d_dt
        }
        
        # --- B. Constraints ---
        if (!is.null(g_dat$geneticConstraint) && length(g_dat$geneticConstraint) > 0) {
          c_dt <- data.table::as.data.table(g_dat$geneticConstraint)
          c_dt[, human_ensembl_id := h_id]
          results_constr[[length(results_constr) + 1]] <- c_dt
        }
        
        # --- C. Tractability ---
        if (!is.null(g_dat$tractability) && length(g_dat$tractability) > 0) {
          t_dt <- data.table::as.data.table(g_dat$tractability)
          t_dt[, human_ensembl_id := h_id]
          results_tract[[length(results_tract) + 1]] <- t_dt
        }
      }
    }, error = function(e) warning("Batch ", i, " failed: ", e$message))
    
    Sys.sleep(0.5)
  }
  
  # 6. Save Results
  save_ot_table <- function(list_data, suffix) {
    if (length(list_data) > 0) {
      dt <- data.table::rbindlist(list_data, fill = TRUE, use.names = TRUE)
      
      # Merge back original symbols (e.g. Mouse Symbols)
      dt <- merge(dt, id_map, by = "human_ensembl_id", all.x = TRUE)
      
      # Clean up: remove disease column if it persists (rare edge case)
      if ("disease" %in% names(dt)) dt[, disease := NULL]
      
      # Reorder to put symbol first
      if ("gene_symbol" %in% names(dt)) {
        col_order <- c("gene_symbol", setdiff(names(dt), "gene_symbol"))
        data.table::setcolorder(dt, col_order)
      }
      
      addRatTable(
        data = dt,
        table_name = paste0("opentargets_", suffix),
        link_type = "gene",
        link_by = "gene_symbol",
        project_dir = project_dir
      )
    }
  }
  
  save_ot_table(results_assoc, "diseases")
  save_ot_table(results_constr, "constraints")
  save_ot_table(results_tract, "tractability")
  
  invisible(list(diseases = results_assoc, constraints = results_constr, tractability = results_tract))
}