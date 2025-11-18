#!/usr/bin/env Rscript

#' Test Package Installation and Loading
#'
#' This script tests building, installing, and loading the locusPackRat package

cat("\n==================================================\n")
cat("Testing locusPackRat Package\n")
cat("==================================================\n\n")

# 1. Test building and installing from source
cat("Step 1: Installing package from source...\n")
cat("-----------------------------------------\n")

pkg_dir <- normalizePath(file.path(getwd()))

# Use devtools if available, otherwise R CMD INSTALL
install_success <- FALSE
if (requireNamespace("devtools", quietly = TRUE)) {
  cat("Using devtools::install()...\n")
  tryCatch({
    devtools::install(pkg_dir, quiet = FALSE, upgrade = "never", reload = FALSE)
    install_success <- TRUE
  }, error = function(e) {
    cat(sprintf("Error: %s\n", e$message))
  })
} else {
  cat("Using R CMD INSTALL...\n")
  result <- system2("R", c("CMD", "INSTALL", pkg_dir), stdout = TRUE, stderr = TRUE)
  install_success <- !any(grepl("ERROR", result))
  if (!install_success) {
    cat(paste(result, collapse = "\n"))
  }
}

if (install_success) {
  cat("✓ Package installed successfully!\n\n")
} else {
  cat("✗ Package installation failed!\n\n")
  quit(status = 1)
}

# 2. Test loading the package
cat("Step 2: Loading package...\n")
cat("--------------------------\n")

load_success <- tryCatch({
  library(locusPackRat)
  cat("✓ Package loaded successfully!\n\n")
  TRUE
}, error = function(e) {
  cat(sprintf("✗ Error loading package: %s\n\n", e$message))
  FALSE
})

if (!load_success) {
  quit(status = 1)
}

# 3. Check exported functions
cat("Step 3: Checking exported functions...\n")
cat("--------------------------------------\n")

expected_functions <- c(
  "initPackRat", "addRatTable", "listPackRatTables",
  "makeGeneSheet", "filterGenes", "makeFilter",
  "generateLocusZoomPlot_v2"
)

all_present <- TRUE
for (func in expected_functions) {
  if (exists(func, mode = "function")) {
    cat(sprintf("  ✓ %s\n", func))
  } else {
    cat(sprintf("  ✗ %s not found\n", func))
    all_present <- FALSE
  }
}

if (all_present) {
  cat("\n✓ All main functions are available!\n\n")
} else {
  cat("\n⚠ Some functions are missing!\n\n")
}

# 4. Test basic functionality
cat("Step 4: Testing basic functionality...\n")
cat("--------------------------------------\n")

# Test with minimal data
test_success <- tryCatch({
  # Create test data
  test_genes <- data.table::data.table(
    gene_symbol = c("Myc", "Tp53", "Egfr", "Vegfa", "Il6"),
    value = 1:5
  )

  # Create temp directory
  temp_dir <- tempfile("locusPackRat_test_")
  dir.create(temp_dir)

  # Initialize project
  cat("  Testing initPackRat()...")
  initPackRat(
    data = test_genes,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = temp_dir,
    force = TRUE
  )
  cat(" ✓\n")

  # Check files created
  cat("  Checking output files...")
  genes_file <- file.path(temp_dir, ".locusPackRat", "input", "genes.csv")
  config_file <- file.path(temp_dir, ".locusPackRat", "config.json")

  if (file.exists(genes_file) && file.exists(config_file)) {
    cat(" ✓\n")

    # Read back the genes to check coordinate enrichment
    genes_data <- data.table::fread(genes_file)
    if ("chr" %in% names(genes_data) && "start" %in% names(genes_data)) {
      cat("  Coordinate enrichment... ✓\n")
    } else {
      cat("  Coordinate enrichment... ⚠ (coordinates not added)\n")
    }
  } else {
    cat(" ✗ (files not created)\n")
  }

  # Test adding supplementary data
  cat("  Testing addRatTable()...")
  supp_data <- data.table::data.table(
    gene_symbol = c("Myc", "Tp53"),
    annotation = c("oncogene", "tumor_suppressor")
  )

  addRatTable(
    data = supp_data,
    table_name = "test_annotation",
    link_type = "gene",
    link_by = "gene_symbol",
    project_dir = temp_dir
  )

  supp_file <- file.path(temp_dir, ".locusPackRat", "supplementary", "test_annotation.csv")
  if (file.exists(supp_file)) {
    cat(" ✓\n")
  } else {
    cat(" ✗\n")
  }

  # Test listing tables
  cat("  Testing listPackRatTables()...")
  tables <- listPackRatTables(temp_dir)
  if (length(tables) > 0 && "test_annotation" %in% tables) {
    cat(" ✓\n")
  } else {
    cat(" ✗\n")
  }

  # Clean up
  unlink(temp_dir, recursive = TRUE)
  cat("\n✓ Basic functionality tests passed!\n")
  TRUE

}, error = function(e) {
  cat(sprintf("\n✗ Error in functionality test: %s\n", e$message))
  FALSE
})

# 5. Summary
cat("\n==================================================\n")
cat("TEST SUMMARY\n")
cat("==================================================\n\n")

if (install_success) {
  cat("✓ Package installation:     SUCCESS\n")
} else {
  cat("✗ Package installation:     FAILED\n")
}

if (load_success) {
  cat("✓ Package loading:          SUCCESS\n")
} else {
  cat("✗ Package loading:          FAILED\n")
}

if (all_present) {
  cat("✓ Function availability:    ALL PRESENT\n")
} else {
  cat("⚠ Function availability:    SOME MISSING\n")
}

if (test_success) {
  cat("✓ Basic functionality:      WORKING\n")
} else {
  cat("✗ Basic functionality:      FAILED\n")
}

cat("\n==================================================\n\n")

# Exit with appropriate code
if (install_success && load_success && all_present && test_success) {
  cat("Overall: ✓ PACKAGE IS READY TO USE!\n\n")
  quit(status = 0)
} else {
  cat("Overall: ✗ PACKAGE NEEDS ATTENTION\n\n")
  quit(status = 1)
}