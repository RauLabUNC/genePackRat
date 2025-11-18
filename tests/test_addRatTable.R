#!/usr/bin/env Rscript

# Test script for addRatTable() function
# Tests adding supplementary tables to existing projects

# Set library path for conda environment
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (conda_prefix != "") {
  conda_lib_path <- file.path(conda_prefix, "lib", "R", "library")
  .libPaths(c(conda_lib_path, .libPaths()))
  cat("Using conda R library path:", conda_lib_path, "\n\n")
}

# Load required libraries
library(data.table)
library(jsonlite)

# Source the functions
source("R/mergeLocusInfo.R")

# Set up test directories
test_base_dir <- file.path(tempdir(), "addRatTable_tests")
dir.create(test_base_dir, showWarnings = FALSE, recursive = TRUE)
cat("Test base directory:", test_base_dir, "\n\n")

# ========================================
# SETUP: Initialize test projects
# ========================================
cat("========================================\n")
cat("SETUP: Creating test projects\n")
cat("========================================\n\n")

# Load test data
test_genes <- fread("tests/test_data/test_gene_list.csv")
test_regions <- fread("tests/test_data/test_regions_chr5.csv")
expression_data <- fread("tests/test_data/test_expression_data.csv")
phenotype_data <- fread("tests/test_data/test_phenotype_data.csv")

# Create gene mode project
gene_project_dir <- file.path(test_base_dir, "gene_project")
dir.create(gene_project_dir, showWarnings = FALSE)

cat("Creating gene mode project...\n")
initPackRat(
  data = test_genes,
  mode = "gene",
  species = "mouse",
  genome = "mm39",
  project_dir = gene_project_dir,
  force = TRUE
)
cat("✓ Gene mode project created\n\n")

# Create region mode project
region_project_dir <- file.path(test_base_dir, "region_project")
dir.create(region_project_dir, showWarnings = FALSE)

cat("Creating region mode project...\n")
initPackRat(
  data = test_regions,
  mode = "region",
  species = "mouse",
  genome = "mm39",
  project_dir = region_project_dir,
  force = TRUE
)
cat("✓ Region mode project created\n\n")

# ========================================
# TEST 1: Add expression data to gene project
# ========================================
cat("========================================\n")
cat("TEST 1: Add Expression Data to Gene Project\n")
cat("========================================\n\n")

cat("Expression data has", nrow(expression_data), "rows\n")
cat("Columns:", paste(names(expression_data), collapse = ", "), "\n\n")

tryCatch({
  # Add expression data - should auto-detect gene_symbol column
  cat("Adding expression data with auto-detection...\n")

  result <- addRatTable(
    data = expression_data,
    table_name = "expression",
    link_type = "gene",
    link_by = NULL,  # Test auto-detection
    project_dir = gene_project_dir
  )

  cat("✓ Expression data added successfully\n\n")

  # Check if file was created
  supp_file <- file.path(gene_project_dir, ".locusPackRat/supplementary/expression.csv")
  if (file.exists(supp_file)) {
    linked_data <- fread(supp_file)
    cat("  ✓ Supplementary file created\n")
    cat("    - Rows:", nrow(linked_data), "\n")
    cat("    - Columns:", paste(names(linked_data), collapse = ", "), "\n")

    # Check if expression columns are present
    if (all(c("expression_ctrl", "expression_treated", "log2fc") %in% names(linked_data))) {
      cat("    - ✓ Expression columns preserved\n")
    }
  } else {
    cat("  ✗ Supplementary file not found\n")
  }

  # Check config update
  config_file <- file.path(gene_project_dir, ".locusPackRat/config.json")
  config <- read_json(config_file)

  if (!is.null(config$supplementary_tables) && "expression" %in% names(config$supplementary_tables)) {
    cat("  ✓ Config updated with table info\n")
    exp_info <- config$supplementary_tables$expression
    cat("    - Link type:", exp_info$link_type, "\n")
    cat("    - Link by:", exp_info$link_by, "\n")
    cat("    - N rows:", exp_info$n_rows, "\n")
    cat("    - N cols:", exp_info$n_cols, "\n")
  }

}, error = function(e) {
  cat("✗ ERROR:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 2: Add phenotype data with explicit link
# ========================================
cat("========================================\n")
cat("TEST 2: Add Phenotype Data with Explicit Link\n")
cat("========================================\n\n")

cat("Phenotype data has", nrow(phenotype_data), "rows\n")
cat("Columns:", paste(names(phenotype_data), collapse = ", "), "\n\n")

tryCatch({
  # Add phenotype data with explicit link_by
  cat("Adding phenotype data with explicit link_by='gene_symbol'...\n")

  result <- addRatTable(
    data = phenotype_data,
    table_name = "phenotypes",
    link_type = "gene",
    link_by = "gene_symbol",  # Explicit linking column
    project_dir = gene_project_dir
  )

  cat("✓ Phenotype data added successfully\n\n")

  # Check if file was created
  supp_file <- file.path(gene_project_dir, ".locusPackRat/supplementary/phenotypes.csv")
  if (file.exists(supp_file)) {
    linked_data <- fread(supp_file)
    cat("  ✓ Supplementary file created\n")
    cat("    - Rows:", nrow(linked_data), "\n")
    cat("    - Columns:", paste(names(linked_data), collapse = ", "), "\n")

    # Check for phenotype columns
    if ("phenotype" %in% names(linked_data) && "p_value" %in% names(linked_data)) {
      cat("    - ✓ Phenotype columns preserved\n")
    }

    # Check for unique phenotypes
    unique_phenotypes <- unique(linked_data$phenotype)
    cat("    - Unique phenotypes:", length(unique_phenotypes), "\n")
    if (length(unique_phenotypes) <= 5) {
      cat("      ", paste(unique_phenotypes, collapse = ", "), "\n")
    }
  } else {
    cat("  ✗ Supplementary file not found\n")
  }

}, error = function(e) {
  cat("✗ ERROR:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 3: Add region-based data
# ========================================
cat("========================================\n")
cat("TEST 3: Add Region-Based Supplementary Data\n")
cat("========================================\n\n")

# Create region-based supplementary data
region_stats <- data.table(
  chr = 5,
  start = c(28340000, 28350000),
  end = c(28360000, 28370000),
  snp_density = c(0.02, 0.03),
  conservation_score = c(0.8, 0.9),
  annotation = c("intergenic", "intronic")
)

cat("Region stats data:\n")
print(region_stats)
cat("\n")

tryCatch({
  cat("Adding region stats to region project...\n")

  result <- addRatTable(
    data = region_stats,
    table_name = "region_stats",
    link_type = "region",
    link_by = NULL,  # Should use coordinate overlap
    project_dir = region_project_dir
  )

  cat("✓ Region stats added successfully\n\n")

  # Check if file was created
  supp_file <- file.path(region_project_dir, ".locusPackRat/supplementary/region_stats.csv")
  if (file.exists(supp_file)) {
    linked_data <- fread(supp_file)
    cat("  ✓ Supplementary file created\n")
    cat("    - Rows:", nrow(linked_data), "\n")
    cat("    - Columns:", paste(names(linked_data), collapse = ", "), "\n")

    if ("snp_density" %in% names(linked_data) && "conservation_score" %in% names(linked_data)) {
      cat("    - ✓ Region statistics preserved\n")
    }
  } else {
    cat("  ✗ Supplementary file not found\n")
  }

}, error = function(e) {
  cat("✗ ERROR:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 4: Test cross-mode linking (genes to region project)
# ========================================
cat("========================================\n")
cat("TEST 4: Cross-Mode Linking\n")
cat("========================================\n\n")

# Since regions don't have genes yet (placeholders), this should fail appropriately
cat("Testing adding gene data to region-mode project...\n")

# Create minimal gene data
gene_data <- data.table(
  gene_symbol = c("Snx18", "Plcb4"),
  gene_function = c("vesicle trafficking", "signal transduction")
)

tryCatch({
  result <- addRatTable(
    data = gene_data,
    table_name = "gene_functions",
    link_type = "gene",
    link_by = "gene_symbol",
    project_dir = region_project_dir
  )

  cat("✓ Gene data added to region project\n")

}, error = function(e) {
  cat("✓ Expected error caught:", e$message, "\n")
  cat("  (This is expected since regions don't have gene mappings yet)\n")
})

cat("\n")

# ========================================
# TEST 5: Error handling
# ========================================
cat("========================================\n")
cat("TEST 5: Error Handling\n")
cat("========================================\n\n")

# Test 5a: Non-existent project directory
cat("Test 5a: Non-existent project directory\n")
tryCatch({
  addRatTable(
    data = expression_data,
    table_name = "test",
    link_type = "gene",
    project_dir = file.path(test_base_dir, "nonexistent")
  )
  cat("  ✗ Should have thrown an error\n")
}, error = function(e) {
  cat("  ✓ Correctly caught error:", e$message, "\n")
})

cat("\n")

# Test 5b: Missing table_name
cat("Test 5b: Missing table_name\n")
tryCatch({
  addRatTable(
    data = expression_data,
    link_type = "gene",
    project_dir = gene_project_dir
  )
  cat("  ✗ Should have thrown an error\n")
}, error = function(e) {
  cat("  ✓ Correctly caught error:", e$message, "\n")
})

cat("\n")

# Test 5c: Non-matching link column
cat("Test 5c: Non-matching link column\n")
bad_data <- data.table(
  wrong_column = c("A", "B"),
  value = c(1, 2)
)

tryCatch({
  addRatTable(
    data = bad_data,
    table_name = "bad_data",
    link_type = "gene",
    link_by = "nonexistent_column",
    project_dir = gene_project_dir
  )
  cat("  ✗ Should have thrown an error\n")
}, error = function(e) {
  cat("  ✓ Correctly caught error:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 6: Multiple tables and merging
# ========================================
cat("========================================\n")
cat("TEST 6: Adding Multiple Tables\n")
cat("========================================\n\n")

# Create another supplementary table
metadata_data <- data.table(
  gene_symbol = test_genes$gene_symbol[1:30],
  pathway = sample(c("metabolism", "signaling", "transport"), 30, replace = TRUE),
  essentiality = sample(c("essential", "non-essential"), 30, replace = TRUE)
)

cat("Adding metadata table...\n")
tryCatch({
  result <- addRatTable(
    data = metadata_data,
    table_name = "gene_metadata",
    link_type = "gene",
    project_dir = gene_project_dir
  )

  cat("✓ Metadata added successfully\n\n")

  # Check all supplementary tables
  supp_dir <- file.path(gene_project_dir, ".locusPackRat/supplementary")
  if (dir.exists(supp_dir)) {
    supp_files <- list.files(supp_dir, pattern = "\\.csv$")
    cat("  All supplementary tables:\n")
    for (f in supp_files) {
      file_size <- file.info(file.path(supp_dir, f))$size
      cat(sprintf("    - %s (%.1f KB)\n", f, file_size/1024))
    }
  }

  # Check config for all tables
  config <- read_json(config_file)
  if (!is.null(config$supplementary_tables)) {
    cat("\n  Tables in config:\n")
    for (tbl_name in names(config$supplementary_tables)) {
      tbl_info <- config$supplementary_tables[[tbl_name]]
      cat(sprintf("    - %s: %d rows, %d cols, linked by %s\n",
                  tbl_name, tbl_info$n_rows, tbl_info$n_cols, tbl_info$link_by))
    }
  }

}, error = function(e) {
  cat("✗ ERROR:", e$message, "\n")
})

cat("\n")

# ========================================
# SUMMARY
# ========================================
cat("========================================\n")
cat("TEST SUMMARY\n")
cat("========================================\n\n")

cat("Test directories:\n")
cat("  - Gene project:", gene_project_dir, "\n")
cat("  - Region project:", region_project_dir, "\n\n")

# Summary of gene project
gene_supp_dir <- file.path(gene_project_dir, ".locusPackRat/supplementary")
if (dir.exists(gene_supp_dir)) {
  gene_tables <- list.files(gene_supp_dir, pattern = "\\.csv$")
  cat("Gene project supplementary tables:", length(gene_tables), "\n")
  for (tbl in gene_tables) {
    cat("  -", tbl, "\n")
  }
}

# Summary of region project
region_supp_dir <- file.path(region_project_dir, ".locusPackRat/supplementary")
if (dir.exists(region_supp_dir)) {
  region_tables <- list.files(region_supp_dir, pattern = "\\.csv$")
  cat("\nRegion project supplementary tables:", length(region_tables), "\n")
  for (tbl in region_tables) {
    cat("  -", tbl, "\n")
  }
}

cat("\nAll tests completed!\n")
cat("To clean up test files, remove:", test_base_dir, "\n")