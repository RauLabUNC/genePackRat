#!/usr/bin/env Rscript

#' Test Package Installation and Dependencies
#'
#' This script tests installing the locusPackRat package in a fresh environment
#' and verifies all dependencies are properly handled.

cat("\n==================================================\n")
cat("Testing locusPackRat Package Installation\n")
cat("==================================================\n\n")

# Function to check if package is installed
check_package <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) "unknown"
    )
    return(list(installed = TRUE, version = version))
  } else {
    return(list(installed = FALSE, version = NA))
  }
}

# List all expected dependencies from DESCRIPTION
imports <- c("data.table", "jsonlite", "openxlsx", "dplyr", "tidyr",
             "plotgardener", "RColorBrewer")
suggests <- c("BiocManager", "GenomicRanges", "GenomeInfoDb", "BiocGenerics",
              "testthat", "knitr", "rmarkdown", "devtools", "usethis",
              "TxDb.Mmusculus.UCSC.mm39.knownGene", "org.Mm.eg.db")

# Check current status of dependencies
cat("Checking current dependency status...\n")
cat("------------------------------------\n\n")

cat("REQUIRED PACKAGES (Imports):\n")
for (pkg in imports) {
  status <- check_package(pkg)
  if (status$installed) {
    cat(sprintf("  ✓ %s (v%s)\n", pkg, status$version))
  } else {
    cat(sprintf("  ✗ %s (not installed)\n", pkg))
  }
}

cat("\nOPTIONAL PACKAGES (Suggests):\n")
for (pkg in suggests) {
  status <- check_package(pkg)
  if (status$installed) {
    cat(sprintf("  ✓ %s (v%s)\n", pkg, status$version))
  } else {
    cat(sprintf("  - %s (not installed)\n", pkg))
  }
}

# Try to install missing dependencies
cat("\n==================================================\n")
cat("Installing missing dependencies...\n")
cat("==================================================\n\n")

# Check for BiocManager (needed for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
}

# Install missing required packages
missing_imports <- imports[!sapply(imports, function(x) requireNamespace(x, quietly = TRUE))]
if (length(missing_imports) > 0) {
  cat("Installing missing required packages:\n")
  for (pkg in missing_imports) {
    cat(sprintf("  Installing %s...\n", pkg))

    # Special handling for plotgardener (Bioconductor package)
    if (pkg == "plotgardener") {
      tryCatch(
        BiocManager::install(pkg, ask = FALSE, quiet = TRUE),
        error = function(e) cat(sprintf("    Error: %s\n", e$message))
      )
    } else {
      tryCatch(
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE),
        error = function(e) cat(sprintf("    Error: %s\n", e$message))
      )
    }
  }
}

# Test package installation
cat("\n==================================================\n")
cat("Testing package build and installation...\n")
cat("==================================================\n\n")

# Get the package directory
pkg_dir <- dirname(dirname(sys.frame(1)$ofile))
if (is.null(pkg_dir)) {
  pkg_dir <- getwd()
}

cat(sprintf("Package directory: %s\n\n", pkg_dir))

# Try to build and install the package
cat("Building package...\n")
build_result <- tryCatch({
  # Use devtools to build and install
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::install(pkg_dir, quiet = FALSE, upgrade = "never")
    TRUE
  } else {
    # Fallback to R CMD
    system2("R", c("CMD", "build", pkg_dir), stdout = TRUE, stderr = TRUE)
    pkg_file <- list.files(pattern = "^locusPackRat.*\\.tar\\.gz$")[1]
    if (!is.na(pkg_file)) {
      install.packages(pkg_file, repos = NULL, type = "source")
      TRUE
    } else {
      FALSE
    }
  }
}, error = function(e) {
  cat(sprintf("Build error: %s\n", e$message))
  FALSE
})

# Test loading the package
cat("\n==================================================\n")
cat("Testing package loading...\n")
cat("==================================================\n\n")

if (build_result) {
  load_test <- tryCatch({
    library(locusPackRat)
    cat("✓ Package loaded successfully!\n\n")

    # Check if main functions are available
    cat("Checking main functions:\n")
    main_functions <- c("initPackRat", "addRatTable", "listPackRatTables",
                       "generateGeneSheet", "filterGenes", "makeFilter")

    for (func in main_functions) {
      if (exists(func)) {
        cat(sprintf("  ✓ %s\n", func))
      } else {
        cat(sprintf("  ✗ %s not found\n", func))
      }
    }

    TRUE
  }, error = function(e) {
    cat(sprintf("✗ Error loading package: %s\n", e$message))
    FALSE
  })
} else {
  cat("✗ Package build failed, cannot test loading\n")
  load_test <- FALSE
}

# Test basic functionality
cat("\n==================================================\n")
cat("Testing basic functionality...\n")
cat("==================================================\n\n")

if (load_test) {
  # Test creating a simple project
  cat("Testing initPackRat with sample data...\n")

  test_result <- tryCatch({
    # Create test data
    test_genes <- data.table::data.table(
      gene_symbol = c("Myc", "Tp53", "Egfr"),
      value = c(1, 2, 3)
    )

    # Create temporary directory for test
    temp_dir <- tempfile("locusPackRat_test")
    dir.create(temp_dir)

    # Initialize project
    initPackRat(
      data = test_genes,
      mode = "gene",
      species = "mouse",
      genome = "mm39",
      project_dir = temp_dir,
      force = TRUE
    )

    cat("✓ initPackRat executed successfully\n")

    # Check if files were created
    expected_files <- c(
      file.path(temp_dir, ".locusPackRat", "input", "genes.csv"),
      file.path(temp_dir, ".locusPackRat", "config.json")
    )

    for (f in expected_files) {
      if (file.exists(f)) {
        cat(sprintf("  ✓ Created: %s\n", basename(f)))
      } else {
        cat(sprintf("  ✗ Missing: %s\n", basename(f)))
      }
    }

    # Clean up
    unlink(temp_dir, recursive = TRUE)
    TRUE

  }, error = function(e) {
    cat(sprintf("✗ Error in basic functionality test: %s\n", e$message))
    FALSE
  })
}

# Final summary
cat("\n==================================================\n")
cat("INSTALLATION TEST SUMMARY\n")
cat("==================================================\n\n")

# Re-check all dependencies after installation attempts
cat("Final dependency check:\n\n")
all_installed <- TRUE

cat("REQUIRED (must have):\n")
for (pkg in imports) {
  status <- check_package(pkg)
  if (status$installed) {
    cat(sprintf("  ✓ %s\n", pkg))
  } else {
    cat(sprintf("  ✗ %s - MISSING (required!)\n", pkg))
    all_installed <- FALSE
  }
}

if (all_installed) {
  cat("\n✓ All required dependencies are installed!\n")
} else {
  cat("\n✗ Some required dependencies are missing!\n")
  cat("  Please install missing packages before using locusPackRat.\n")
}

if (exists("load_test") && load_test) {
  cat("\n✓ Package loads successfully!\n")
} else {
  cat("\n✗ Package failed to load properly.\n")
}

if (exists("test_result") && test_result) {
  cat("✓ Basic functionality test passed!\n")
} else {
  cat("✗ Basic functionality test failed.\n")
}

cat("\n==================================================\n\n")