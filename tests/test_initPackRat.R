#!/usr/bin/env Rscript

# Test script for initPackRat() function
# Tests both gene mode and region mode initialization

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
test_base_dir <- file.path(tempdir(), "initPackRat_tests")
dir.create(test_base_dir, showWarnings = FALSE, recursive = TRUE)
cat("Test base directory:", test_base_dir, "\n\n")

# ========================================
# TEST 1: Gene Mode Initialization
# ========================================
cat("========================================\n")
cat("TEST 1: Gene Mode Initialization\n")
cat("========================================\n\n")

# Load test gene data
test_genes <- fread("tests/test_data/test_gene_list.csv")
cat("Loaded", nrow(test_genes), "test genes\n")
cat("Columns:", paste(names(test_genes), collapse = ", "), "\n\n")

# Create test directory for gene mode
gene_test_dir <- file.path(test_base_dir, "gene_mode")
dir.create(gene_test_dir, showWarnings = FALSE)

# Test initialization
cat("Initializing project in gene mode...\n")
tryCatch({
  result <- initPackRat(
    data = test_genes,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = gene_test_dir,
    force = TRUE
  )

  cat("✓ initPackRat() completed successfully\n\n")

  # Verify created structure
  packrat_dir <- file.path(gene_test_dir, ".locusPackRat")

  # Check directories
  expected_dirs <- c("input", "supplementary", "cache", "output")
  cat("Checking directory structure:\n")
  for (d in expected_dirs) {
    dir_path <- file.path(packrat_dir, d)
    if (dir.exists(dir_path)) {
      cat("  ✓", d, "directory created\n")
    } else {
      cat("  ✗", d, "directory missing\n")
    }
  }

  # Check files
  cat("\nChecking created files:\n")

  # Config file
  config_file <- file.path(packrat_dir, "config.json")
  if (file.exists(config_file)) {
    cat("  ✓ config.json created\n")
    config <- read_json(config_file)
    cat("    - Mode:", config$mode, "\n")
    cat("    - Species:", config$species, "\n")
    cat("    - Genome:", config$genome, "\n")
    cat("    - N entries:", config$n_entries, "\n")
    cat("    - Input columns:", paste(config$input_columns, collapse = ", "), "\n")
  } else {
    cat("  ✗ config.json missing\n")
  }

  # Genes file
  genes_file <- file.path(packrat_dir, "input/genes.csv")
  if (file.exists(genes_file)) {
    genes_data <- fread(genes_file)
    cat("  ✓ genes.csv created (", nrow(genes_data), "rows )\n")
    cat("    - Columns:", paste(names(genes_data), collapse = ", "), "\n")

    # Check for preserved custom columns
    if ("custom_value" %in% names(genes_data)) {
      cat("    - ✓ Custom column 'custom_value' preserved\n")
    }
    if ("category" %in% names(genes_data)) {
      cat("    - ✓ Custom column 'category' preserved\n")
    }
  } else {
    cat("  ✗ genes.csv missing\n")
  }

  # Orthology file
  ortho_file <- file.path(packrat_dir, "input/orthology.csv")
  if (file.exists(ortho_file)) {
    ortho_data <- fread(ortho_file)
    cat("  ✓ orthology.csv created (", nrow(ortho_data), "rows )\n")
    cat("    - Columns:", paste(names(ortho_data), collapse = ", "), "\n")
  } else {
    cat("  ✗ orthology.csv missing\n")
  }

}, error = function(e) {
  cat("✗ ERROR:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 2: Region Mode Initialization
# ========================================
cat("========================================\n")
cat("TEST 2: Region Mode Initialization\n")
cat("========================================\n\n")

# Load test region data
test_regions <- fread("tests/test_data/test_regions_chr5.csv")
cat("Loaded", nrow(test_regions), "test regions from chromosome 5\n")
cat("Columns:", paste(names(test_regions), collapse = ", "), "\n\n")

# Create test directory for region mode
region_test_dir <- file.path(test_base_dir, "region_mode")
dir.create(region_test_dir, showWarnings = FALSE)

# Test initialization
cat("Initializing project in region mode...\n")
tryCatch({
  result <- initPackRat(
    data = test_regions,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = region_test_dir,
    force = TRUE
  )

  cat("✓ initPackRat() completed successfully\n\n")

  # Verify created structure
  packrat_dir <- file.path(region_test_dir, ".locusPackRat")

  # Check directories
  expected_dirs <- c("input", "supplementary", "cache", "output")
  cat("Checking directory structure:\n")
  for (d in expected_dirs) {
    dir_path <- file.path(packrat_dir, d)
    if (dir.exists(dir_path)) {
      cat("  ✓", d, "directory created\n")
    } else {
      cat("  ✗", d, "directory missing\n")
    }
  }

  # Check files
  cat("\nChecking created files:\n")

  # Config file
  config_file <- file.path(packrat_dir, "config.json")
  if (file.exists(config_file)) {
    cat("  ✓ config.json created\n")
    config <- read_json(config_file)
    cat("    - Mode:", config$mode, "\n")
    cat("    - Species:", config$species, "\n")
    cat("    - Genome:", config$genome, "\n")
    cat("    - N entries:", config$n_entries, "\n")
  } else {
    cat("  ✗ config.json missing\n")
  }

  # Regions file
  regions_file <- file.path(packrat_dir, "input/regions.csv")
  if (file.exists(regions_file)) {
    regions_data <- fread(regions_file)
    cat("  ✓ regions.csv created (", nrow(regions_data), "rows )\n")
    cat("    - Columns:", paste(names(regions_data), collapse = ", "), "\n")

    # Check for preserved custom columns
    if ("peak_pos" %in% names(regions_data)) {
      cat("    - ✓ Custom column 'peak_pos' preserved\n")
    }
    if ("max_lod" %in% names(regions_data)) {
      cat("    - ✓ Custom column 'max_lod' preserved\n")
    }
    if ("trait" %in% names(regions_data)) {
      cat("    - ✓ Custom column 'trait' preserved\n")
    }
    if ("drug" %in% names(regions_data)) {
      cat("    - ✓ Custom column 'drug' preserved\n")
    }

    # Check coordinate values
    sample_region <- regions_data[1]
    cat("    - Sample region: chr", sample_region$chr,
        " [", format(sample_region$start, scientific = FALSE),
        "-", format(sample_region$end, scientific = FALSE), "]\n")
  } else {
    cat("  ✗ regions.csv missing\n")
  }

  # Orthology file (should exist even for regions)
  ortho_file <- file.path(packrat_dir, "input/orthology.csv")
  if (file.exists(ortho_file)) {
    ortho_data <- fread(ortho_file)
    cat("  ✓ orthology.csv created (", nrow(ortho_data), "rows )\n")
  } else {
    cat("  ✗ orthology.csv missing (expected for regions mode)\n")
  }

}, error = function(e) {
  cat("✗ ERROR:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 3: Error Handling
# ========================================
cat("========================================\n")
cat("TEST 3: Error Handling\n")
cat("========================================\n\n")

# Test 3a: Missing required columns
cat("Test 3a: Missing required columns for regions\n")
bad_regions <- data.table(
  chromosome = 1,  # Wrong column name
  start_pos = 1000,
  end_pos = 2000
)

error_test_dir <- file.path(test_base_dir, "error_test")
dir.create(error_test_dir, showWarnings = FALSE)

tryCatch({
  initPackRat(
    data = bad_regions,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = error_test_dir,
    force = TRUE
  )
  cat("  ✗ Should have thrown an error for missing columns\n")
}, error = function(e) {
  cat("  ✓ Correctly caught error:", e$message, "\n")
})

cat("\n")

# Test 3b: Incompatible species/genome
cat("Test 3b: Incompatible species/genome combination\n")
tryCatch({
  initPackRat(
    data = test_genes,
    mode = "gene",
    species = "human",
    genome = "mm39",  # Mouse genome with human species
    project_dir = error_test_dir,
    force = TRUE
  )
  cat("  ✗ Should have thrown an error for incompatible species/genome\n")
}, error = function(e) {
  cat("  ✓ Correctly caught error:", e$message, "\n")
})

cat("\n")

# Test 3c: Existing directory without force
cat("Test 3c: Existing directory without force=TRUE\n")
# First create a project
existing_dir <- file.path(test_base_dir, "existing_project")
dir.create(existing_dir, showWarnings = FALSE)
initPackRat(
  data = test_genes[1:5],
  mode = "gene",
  species = "mouse",
  genome = "mm39",
  project_dir = existing_dir,
  force = TRUE
)

# Try to initialize again without force
tryCatch({
  initPackRat(
    data = test_genes[1:5],
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = existing_dir,
    force = FALSE
  )
  cat("  ✗ Should have thrown an error for existing directory\n")
}, error = function(e) {
  cat("  ✓ Correctly caught error:", e$message, "\n")
})

cat("\n")

# ========================================
# TEST 4: Human Gene Mode
# ========================================
cat("========================================\n")
cat("TEST 4: Human Gene Mode\n")
cat("========================================\n\n")

# Create test human genes
human_genes <- data.table(
  gene_symbol = c("BRCA1", "TP53", "EGFR", "MYC", "PTEN"),
  expression = runif(5, 0, 100)
)

human_test_dir <- file.path(test_base_dir, "human_mode")
dir.create(human_test_dir, showWarnings = FALSE)

cat("Testing with human genes...\n")
tryCatch({
  result <- initPackRat(
    data = human_genes,
    mode = "gene",
    species = "human",
    genome = "hg38",
    project_dir = human_test_dir,
    force = TRUE
  )

  cat("✓ Human gene mode initialized successfully\n")

  # Check orthology columns
  ortho_file <- file.path(human_test_dir, ".locusPackRat/input/orthology.csv")
  if (file.exists(ortho_file)) {
    ortho_data <- fread(ortho_file)
    cat("  Orthology columns:", paste(names(ortho_data), collapse = ", "), "\n")
    if ("human_gene_symbol" %in% names(ortho_data) &&
        "mouse_gene_symbol" %in% names(ortho_data)) {
      cat("  ✓ Correct species-specific column names\n")
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

cat("Test directories created in:", test_base_dir, "\n")
cat("Projects initialized:\n")

# List all created projects
projects <- list.dirs(test_base_dir, recursive = FALSE)
for (proj in projects) {
  packrat_path <- file.path(proj, ".locusPackRat")
  if (dir.exists(packrat_path)) {
    config_file <- file.path(packrat_path, "config.json")
    if (file.exists(config_file)) {
      config <- read_json(config_file)
      cat("  -", basename(proj), ":", config$mode, "mode,",
          config$species, config$genome, "\n")
    }
  }
}

cat("\nAll tests completed!\n")
cat("To clean up test files, remove:", test_base_dir, "\n")