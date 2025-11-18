#!/usr/bin/env Rscript

# Test script for mergeLocusInfo.R functions
# Tests initPackRat, addRatTable, and generateGeneSheet

library(data.table)
library(jsonlite)
library(openxlsx)

# Source the functions
source("R/mergeLocusInfo.R")

# Set up test environment
test_dir <- tempdir()
cat("Using test directory:", test_dir, "\n\n")

# ========================================
# STEP 1: Create stand-in genome tables
# ========================================
cat("=== Creating stand-in genome coordinate tables ===\n")

# Load and prepare gene coordinate data
gene_info <- fread("inst/extdata/gene_info.csv")
orthology <- fread("inst/extdata/orthology.csv")

# Create mouse genome coordinate table (mm39)
mouse_coords <- gene_info[!is.na(chromosome_name) & chromosome_name != "", .(
  gene_symbol = mouse_gene_symbol,
  ensembl_id = mouse_ensembl_id,
  chr = chromosome_name,
  start = start_position,
  end = end_position,
  strand = "+"  # Placeholder
)]

# Remove duplicates and NA gene symbols
mouse_coords <- mouse_coords[!is.na(gene_symbol) & gene_symbol != ""]
mouse_coords <- unique(mouse_coords, by = "gene_symbol")

cat("Created mouse coordinate table with", nrow(mouse_coords), "genes\n")

# Create human coordinate table (hg38) - using orthology to get human genes
human_genes <- unique(orthology[!is.na(human_gene_symbol) & human_gene_symbol != "", .(
  gene_symbol = human_gene_symbol,
  ensembl_id = human_ensembl_id
)])

# For testing, create dummy coordinates for human genes
set.seed(123)
human_coords <- human_genes[, .(
  gene_symbol,
  ensembl_id,
  chr = sample(c(1:22, "X", "Y"), .N, replace = TRUE),
  start = sample(1000000:100000000, .N),
  end = NA_integer_,
  strand = "+"
)]
human_coords[, end := start + sample(1000:50000, .N)]

cat("Created human coordinate table with", nrow(human_coords), "genes\n")

# Create combined orthology table
ortho_table <- merge(
  orthology[, .(mouse_ensembl_id, human_ensembl_id, human_gene_symbol)],
  gene_info[, .(mouse_ensembl_id, mouse_gene_symbol)],
  by = "mouse_ensembl_id",
  all.x = TRUE
)

cat("Created orthology table with", nrow(ortho_table), "mappings\n\n")

# ========================================
# TEST 1: initPackRat - Gene Mode
# ========================================
cat("=== TEST 1: initPackRat with Gene Mode ===\n")

# Prepare test gene data
test_genes <- mouse_coords[sample(.N, 50), .(
  gene_symbol = gene_symbol,
  my_custom_value = rnorm(50)
)]

# Test with gene mode
gene_test_dir <- file.path(test_dir, "gene_mode_test")
dir.create(gene_test_dir, showWarnings = FALSE)

tryCatch({
  initPackRat(
    data = test_genes,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = gene_test_dir,
    force = TRUE
  )

  # Check if files were created
  packrat_dir <- file.path(gene_test_dir, ".locusPackRat")
  if (dir.exists(packrat_dir)) {
    cat("✓ .locusPackRat directory created\n")

    # Check for expected files
    expected_files <- c(
      "config.json",
      "input/genes.csv",
      "input/orthology.csv"
    )

    for (f in expected_files) {
      fpath <- file.path(packrat_dir, f)
      if (file.exists(fpath)) {
        cat("✓ File created:", f, "\n")
      } else {
        cat("✗ Missing file:", f, "\n")
      }
    }

    # Read and display config
    config <- read_json(file.path(packrat_dir, "config.json"))
    cat("\nConfig contents:\n")
    cat("  Mode:", config$mode, "\n")
    cat("  Species:", config$species, "\n")
    cat("  Genome:", config$genome, "\n")
    cat("  N entries:", config$n_entries, "\n")
  }

}, error = function(e) {
  cat("ERROR in gene mode test:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 2: initPackRat - Region Mode
# ========================================
cat("=== TEST 2: initPackRat with Region Mode ===\n")

# Load and prepare region data
regions <- fread("inst/extdata/all_significant_regions_summary.csv")

# Filter to chromosome 5 and convert column names
test_regions <- regions[chr == "5", .(
  chr = chr,
  start = upper_pos_lod_drop * 1e6,  # Convert to base pairs
  end = lower_pos_lod_drop * 1e6,
  peak_pos = peak_pos,
  max_lod = max_lod,
  trait = trait,
  drug = drug
)]

# Make sure start < end
test_regions[start > end, `:=`(temp = start, start = end, end = temp)]
test_regions[, temp := NULL]

cat("Prepared", nrow(test_regions), "regions from chromosome 5\n")

# Test with region mode
region_test_dir <- file.path(test_dir, "region_mode_test")
dir.create(region_test_dir, showWarnings = FALSE)

tryCatch({
  initPackRat(
    data = test_regions,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = region_test_dir,
    force = TRUE
  )

  # Check if files were created
  packrat_dir <- file.path(region_test_dir, ".locusPackRat")
  if (dir.exists(packrat_dir)) {
    cat("✓ .locusPackRat directory created\n")

    # Check regions file
    regions_file <- file.path(packrat_dir, "input/regions.csv")
    if (file.exists(regions_file)) {
      cat("✓ Regions file created\n")
      saved_regions <- fread(regions_file)
      cat("  Saved", nrow(saved_regions), "regions\n")
      cat("  Columns:", paste(names(saved_regions), collapse = ", "), "\n")
    }
  }

}, error = function(e) {
  cat("ERROR in region mode test:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 3: addRatTable
# ========================================
cat("=== TEST 3: addRatTable ===\n")

# Create supplementary expression data
expression_data <- data.table(
  gene_symbol = test_genes$gene_symbol[1:20],
  expression_level = runif(20, 0, 100),
  condition = sample(c("control", "treated"), 20, replace = TRUE)
)

tryCatch({
  addRatTable(
    data = expression_data,
    table_name = "expression_test",
    link_type = "gene",
    link_by = "gene_symbol",
    project_dir = gene_test_dir
  )

  # Check if supplementary file was created
  supp_file <- file.path(gene_test_dir, ".locusPackRat/supplementary/expression_test.csv")
  if (file.exists(supp_file)) {
    cat("✓ Supplementary table created\n")
    supp_data <- fread(supp_file)
    cat("  Linked", nrow(supp_data), "rows\n")
    cat("  Columns:", paste(names(supp_data), collapse = ", "), "\n")
  }

}, error = function(e) {
  cat("ERROR in addRatTable test:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 4: generateGeneSheet - CSV
# ========================================
cat("=== TEST 4: generateGeneSheet (CSV) ===\n")

tryCatch({
  # Generate CSV output
  csv_output <- file.path(gene_test_dir, ".locusPackRat/output/test_sheet.csv")

  result <- generateGeneSheet(
    filter_expr = NULL,  # No filter for now
    highlight_genes = test_genes$gene_symbol[1:5],
    format = "csv",
    output_file = csv_output,
    include_supplementary = TRUE,
    project_dir = gene_test_dir
  )

  if (file.exists(csv_output)) {
    cat("✓ CSV file created:", csv_output, "\n")
    csv_data <- fread(csv_output)
    cat("  Rows:", nrow(csv_data), "\n")
    cat("  Columns:", ncol(csv_data), "\n")
  }

}, error = function(e) {
  cat("ERROR in generateGeneSheet CSV test:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 5: generateGeneSheet - Excel
# ========================================
cat("=== TEST 5: generateGeneSheet (Excel) ===\n")

tryCatch({
  # Generate Excel output with multiple tabs
  excel_output <- file.path(gene_test_dir, ".locusPackRat/output/test_sheet.xlsx")

  result <- generateGeneSheet(
    filter_expr = NULL,
    highlight_genes = test_genes$gene_symbol[1:5],
    highlight_color = "#FFCCCC",
    format = "excel",
    split_by = "criteria",
    split_criteria = list(
      "High_Value" = "my_custom_value > 0",
      "Low_Value" = "my_custom_value <= 0"
    ),
    output_file = excel_output,
    include_supplementary = TRUE,
    project_dir = gene_test_dir
  )

  if (file.exists(excel_output)) {
    cat("✓ Excel file created:", excel_output, "\n")

    # Check Excel contents
    wb <- loadWorkbook(excel_output)
    sheets <- names(wb)
    cat("  Sheets:", paste(sheets, collapse = ", "), "\n")

    for (sheet in sheets) {
      sheet_data <- read.xlsx(excel_output, sheet = sheet)
      cat("    -", sheet, ":", nrow(sheet_data), "rows\n")
    }
  }

}, error = function(e) {
  cat("ERROR in generateGeneSheet Excel test:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 6: Region Mode Integration
# ========================================
cat("=== TEST 6: Region Mode with generateGeneSheet ===\n")

tryCatch({
  # Test generateGeneSheet with region mode project
  result <- generateGeneSheet(
    summary_level = "locus",
    format = "csv",
    output_file = file.path(region_test_dir, ".locusPackRat/output/region_summary.csv"),
    project_dir = region_test_dir
  )

  output_file <- file.path(region_test_dir, ".locusPackRat/output/region_summary.csv")
  if (file.exists(output_file)) {
    cat("✓ Region summary created\n")
    summary_data <- fread(output_file)
    cat("  Rows:", nrow(summary_data), "\n")
    cat("  Columns:", paste(names(summary_data), collapse = ", "), "\n")
  }

}, error = function(e) {
  cat("ERROR in region mode generateGeneSheet test:", e$message, "\n")
})

cat("\n")

# ========================================
# SUMMARY
# ========================================
cat("=== TEST SUMMARY ===\n")
cat("Test directory:", test_dir, "\n")
cat("Gene mode project:", file.path(gene_test_dir, ".locusPackRat"), "\n")
cat("Region mode project:", file.path(region_test_dir, ".locusPackRat"), "\n")
cat("\nTo explore results, navigate to the test directories listed above.\n")
cat("To clean up, remove the temp directory:", test_dir, "\n")