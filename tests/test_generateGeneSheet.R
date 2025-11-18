#!/usr/bin/env Rscript

# Test script for generateGeneSheet() function
# Tests generating CSV and Excel outputs with various options

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

# Check for openxlsx
has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
if (has_openxlsx) {
  library(openxlsx)
  cat("openxlsx available - Excel tests will be run\n\n")
} else {
  cat("openxlsx not available - Excel tests will be skipped\n\n")
}

# Source the functions
source("R/mergeLocusInfo.R")

# Set up test directories
test_base_dir <- file.path("tests", "output_examples")  # Save in repo for review
dir.create(test_base_dir, showWarnings = FALSE, recursive = TRUE)
cat("Test output directory:", test_base_dir, "\n\n")

# ========================================
# SETUP: Create test project with data
# ========================================
cat("========================================\n")
cat("SETUP: Creating test project with supplementary data\n")
cat("========================================\n\n")

# Load test data
test_genes <- fread("tests/test_data/test_gene_list.csv")
expression_data <- fread("tests/test_data/test_expression_data.csv")
phenotype_data <- fread("tests/test_data/test_phenotype_data.csv")

# Create project directory
project_dir <- file.path(test_base_dir, "test_project")
dir.create(project_dir, showWarnings = FALSE)

# Initialize project
cat("Initializing gene project...\n")
initPackRat(
  data = test_genes,
  mode = "gene",
  species = "mouse",
  genome = "mm39",
  project_dir = project_dir,
  force = TRUE
)

# Add supplementary tables
cat("Adding expression data...\n")
addRatTable(
  data = expression_data,
  table_name = "expression",
  link_type = "gene",
  link_by = "gene_symbol",
  project_dir = project_dir
)

cat("Adding phenotype data...\n")
addRatTable(
  data = phenotype_data,
  table_name = "phenotypes",
  link_type = "gene",
  link_by = "gene_symbol",
  project_dir = project_dir
)

cat("✓ Project setup complete with supplementary data\n\n")

# ========================================
# TEST 1: Basic CSV output
# ========================================
cat("========================================\n")
cat("TEST 1: Basic CSV Output\n")
cat("========================================\n\n")

tryCatch({
  cat("Generating basic CSV gene sheet...\n")

  result <- generateGeneSheet(
    format = "csv",
    include_supplementary = TRUE,
    project_dir = project_dir
  )

  csv_file <- file.path(project_dir, ".locusPackRat/output/gene_sheet.csv")

  if (file.exists(csv_file)) {
    csv_data <- fread(csv_file)
    cat("✓ CSV file created:", csv_file, "\n")
    cat("  - Rows:", nrow(csv_data), "\n")
    cat("  - Columns:", ncol(csv_data), "\n")
    cat("  - Column names:", paste(names(csv_data)[1:min(10, ncol(csv_data))], collapse = ", "))
    if (ncol(csv_data) > 10) cat(", ...")
    cat("\n")

    # Check if supplementary data was included
    if ("expression_ctrl" %in% names(csv_data)) {
      cat("  - ✓ Expression data included\n")
    }
    if ("phenotype" %in% names(csv_data)) {
      cat("  - ✓ Phenotype data included\n")
    }
  } else {
    cat("✗ CSV file not created\n")
  }

}, error = function(e) {
  cat("✗ ERROR:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 2: Excel with highlighting
# ========================================
if (has_openxlsx) {
  cat("========================================\n")
  cat("TEST 2: Excel Output with Highlighting\n")
  cat("========================================\n\n")

  tryCatch({
    cat("Generating Excel with highlighted genes...\n")

    # Get some actual genes from the test data to highlight
    highlight_genes <- c("Snx18", "Plcb4", "Mafa", "Adra2b", "Nckap5")
    cat("Genes to highlight:", paste(highlight_genes, collapse = ", "), "\n")

    result <- generateGeneSheet(
      format = "excel",
      highlight_genes = highlight_genes,
      highlight_color = "#FFCCCC",  # Light red
      include_supplementary = TRUE,
      output_file = file.path(test_base_dir, "example_highlighted.xlsx"),
      project_dir = project_dir
    )

    excel_file <- file.path(test_base_dir, "example_highlighted.xlsx")

    if (file.exists(excel_file)) {
      cat("✓ Excel file created:", excel_file, "\n")

      # Read and check the workbook
      wb <- loadWorkbook(excel_file)
      sheets <- names(wb)
      cat("  - Sheets:", paste(sheets, collapse = ", "), "\n")

      # Check first sheet
      sheet_data <- read.xlsx(excel_file, sheet = 1)
      cat("  - First sheet rows:", nrow(sheet_data), "\n")
      cat("  - First sheet columns:", ncol(sheet_data), "\n")
    } else {
      cat("✗ Excel file not created\n")
    }

  }, error = function(e) {
    cat("✗ ERROR:", e$message, "\n")
  })

  cat("\n")
}

# ========================================
# TEST 3: Filtered output
# ========================================
cat("========================================\n")
cat("TEST 3: Filtered Gene Sheet\n")
cat("========================================\n\n")

tryCatch({
  cat("Generating filtered gene sheet (category == 'A')...\n")

  result <- generateGeneSheet(
    filter_expr = "category == 'A'",
    format = "csv",
    output_file = file.path(test_base_dir, "filtered_genes.csv"),
    include_supplementary = TRUE,
    project_dir = project_dir
  )

  filtered_file <- file.path(test_base_dir, "filtered_genes.csv")

  if (file.exists(filtered_file)) {
    filtered_data <- fread(filtered_file)
    cat("✓ Filtered CSV created:", filtered_file, "\n")
    cat("  - Rows after filter:", nrow(filtered_data), "\n")

    # Check if all are category A
    if ("category" %in% names(filtered_data)) {
      unique_cats <- unique(filtered_data$category)
      cat("  - Categories in filtered data:", paste(unique_cats, collapse = ", "), "\n")
      if (all(unique_cats == "A")) {
        cat("  - ✓ Filter correctly applied\n")
      }
    }
  } else {
    cat("✗ Filtered file not created\n")
  }

}, error = function(e) {
  cat("✗ ERROR:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 4: Excel with multiple tabs
# ========================================
if (has_openxlsx) {
  cat("========================================\n")
  cat("TEST 4: Excel with Multiple Tabs (Split by Criteria)\n")
  cat("========================================\n\n")

  tryCatch({
    cat("Generating Excel with custom tabs based on criteria...\n")

    result <- generateGeneSheet(
      format = "excel",
      split_by = "criteria",
      split_criteria = list(
        "Category_A" = "category == 'A'",
        "Category_B" = "category == 'B'",
        "High_Expression" = "!is.na(expression_ctrl) & expression_ctrl > 50"
      ),
      include_supplementary = TRUE,
      output_file = file.path(test_base_dir, "example_multi_tab.xlsx"),
      project_dir = project_dir
    )

    excel_file <- file.path(test_base_dir, "example_multi_tab.xlsx")

    if (file.exists(excel_file)) {
      cat("✓ Multi-tab Excel created:", excel_file, "\n")

      wb <- loadWorkbook(excel_file)
      sheets <- names(wb)
      cat("  - Sheets created:", paste(sheets, collapse = ", "), "\n")

      # Check each sheet
      for (sheet in sheets) {
        sheet_data <- read.xlsx(excel_file, sheet = sheet)
        cat("    -", sheet, ":", nrow(sheet_data), "rows\n")
      }
    } else {
      cat("✗ Multi-tab Excel not created\n")
    }

  }, error = function(e) {
    cat("✗ ERROR:", e$message, "\n")
  })

  cat("\n")
}

# ========================================
# TEST 5: Highlighting with expression
# ========================================
if (has_openxlsx) {
  cat("========================================\n")
  cat("TEST 5: Excel with Expression-Based Highlighting\n")
  cat("========================================\n\n")

  tryCatch({
    cat("Generating Excel with expression-based highlighting...\n")

    result <- generateGeneSheet(
      format = "excel",
      highlight_expr = "!is.na(p_value) & p_value < 0.05",
      highlight_color = "#CCFFCC",  # Light green
      include_supplementary = TRUE,
      output_file = file.path(test_base_dir, "example_expression_highlight.xlsx"),
      project_dir = project_dir
    )

    excel_file <- file.path(test_base_dir, "example_expression_highlight.xlsx")

    if (file.exists(excel_file)) {
      cat("✓ Expression-highlighted Excel created:", excel_file, "\n")

      sheet_data <- read.xlsx(excel_file, sheet = 1)
      if ("p_value" %in% names(sheet_data)) {
        sig_count <- sum(sheet_data$p_value < 0.05, na.rm = TRUE)
        cat("  - Genes with p_value < 0.05:", sig_count, "(highlighted)\n")
      }
    } else {
      cat("✗ Expression-highlighted Excel not created\n")
    }

  }, error = function(e) {
    cat("✗ ERROR:", e$message, "\n")
  })

  cat("\n")
}

# ========================================
# TEST 6: Summary at locus level
# ========================================
cat("========================================\n")
cat("TEST 6: Summary at Locus Level\n")
cat("========================================\n\n")

tryCatch({
  cat("Generating locus-level summary...\n")

  result <- generateGeneSheet(
    summary_level = "locus",
    format = "csv",
    output_file = file.path(test_base_dir, "locus_summary.csv"),
    include_supplementary = FALSE,  # Simpler without supplementary
    project_dir = project_dir
  )

  summary_file <- file.path(test_base_dir, "locus_summary.csv")

  if (file.exists(summary_file)) {
    summary_data <- fread(summary_file)
    cat("✓ Locus summary created:", summary_file, "\n")
    cat("  - Summary rows:", nrow(summary_data), "\n")
    cat("  - Columns:", paste(names(summary_data), collapse = ", "), "\n")

    if ("n_genes" %in% names(summary_data)) {
      cat("  - Total genes summarized:", sum(summary_data$n_genes), "\n")
    }
  } else {
    cat("✗ Locus summary not created\n")
  }

}, error = function(e) {
  cat("✗ ERROR:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 7: Complete showcase Excel
# ========================================
if (has_openxlsx) {
  cat("========================================\n")
  cat("TEST 7: Complete Showcase Excel\n")
  cat("========================================\n\n")

  tryCatch({
    cat("Generating comprehensive Excel showcase...\n")

    # Combine multiple features
    result <- generateGeneSheet(
      format = "excel",
      highlight_genes = c("Snx18", "Plcb4", "Gapdh"),
      highlight_expr = "!is.na(p_value) & p_value < 0.05",
      highlight_color = "#FFE6CC",  # Light orange
      split_by = "criteria",
      split_criteria = list(
        "Significant" = "!is.na(p_value) & p_value < 0.05",
        "High_Expression" = "!is.na(expression_ctrl) & expression_ctrl > 70",
        "With_Phenotype" = "!is.na(phenotype)"
      ),
      include_supplementary = TRUE,
      output_file = file.path(test_base_dir, "SHOWCASE_complete_example.xlsx"),
      project_dir = project_dir
    )

    excel_file <- file.path(test_base_dir, "SHOWCASE_complete_example.xlsx")

    if (file.exists(excel_file)) {
      cat("✓ SHOWCASE Excel created:", excel_file, "\n")
      cat("  *** This file demonstrates all features ***\n")

      wb <- loadWorkbook(excel_file)
      sheets <- names(wb)
      cat("  - Sheets:", paste(sheets, collapse = ", "), "\n")

      for (sheet in sheets) {
        sheet_data <- read.xlsx(excel_file, sheet = sheet)
        cat("    -", sheet, ":", nrow(sheet_data), "rows,", ncol(sheet_data), "columns\n")
      }

      cat("\n  Features demonstrated:\n")
      cat("    - Multiple tabs based on criteria\n")
      cat("    - Gene highlighting (specific genes + expression)\n")
      cat("    - Supplementary data integration\n")
      cat("    - Formatted headers and striped rows\n")
      cat("    - Auto-filtered tables\n")
    } else {
      cat("✗ Showcase Excel not created\n")
    }

  }, error = function(e) {
    cat("✗ ERROR:", e$message, "\n")
  })

  cat("\n")
}

# ========================================
# SUMMARY
# ========================================
cat("========================================\n")
cat("TEST SUMMARY\n")
cat("========================================\n\n")

cat("Output files created in:", test_base_dir, "\n\n")

# List all created files
output_files <- list.files(test_base_dir, pattern = "\\.(csv|xlsx)$",
                           recursive = FALSE, full.names = FALSE)

if (length(output_files) > 0) {
  cat("Files available for review:\n")
  for (f in output_files) {
    file_path <- file.path(test_base_dir, f)
    file_size <- file.info(file_path)$size
    cat(sprintf("  - %s (%.1f KB)\n", f, file_size/1024))
  }

  cat("\n*** KEY FILE FOR REVIEW ***\n")
  cat("  SHOWCASE_complete_example.xlsx - Demonstrates all Excel features\n")
} else {
  cat("No output files were created\n")
}

cat("\nProject directory:", project_dir, "\n")
cat("Output directory:", test_base_dir, "\n")
cat("\nAll tests completed!\n")