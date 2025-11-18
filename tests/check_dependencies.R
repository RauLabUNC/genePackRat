#!/usr/bin/env Rscript

#' Check Package Dependencies Status
#'
#' Quick script to check which dependencies are installed

cat("\n==================================================\n")
cat("Checking locusPackRat Dependencies\n")
cat("==================================================\n\n")

# Required packages from DESCRIPTION
imports <- c("data.table", "jsonlite", "openxlsx", "dplyr", "tidyr",
             "plotgardener", "RColorBrewer")

cat("REQUIRED PACKAGES (Imports):\n")
cat("--------------------------\n")
for (pkg in imports) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    cat(sprintf("  ✓ %-20s v%s\n", pkg, version))
  } else {
    cat(sprintf("  ✗ %-20s NOT INSTALLED\n", pkg))
  }
}

# Count missing
missing <- imports[!sapply(imports, function(x) requireNamespace(x, quietly = TRUE))]
if (length(missing) > 0) {
  cat("\n⚠ Missing packages:", paste(missing, collapse = ", "), "\n")
  cat("\nTo install missing packages, run:\n")
  if ("plotgardener" %in% missing) {
    cat('  BiocManager::install("plotgardener")\n')
  }
  other_missing <- setdiff(missing, "plotgardener")
  if (length(other_missing) > 0) {
    cat('  install.packages(c("', paste(other_missing, collapse = '", "'), '"))\n', sep = "")
  }
} else {
  cat("\n✓ All required packages are installed!\n")
}

# Check if we're in a conda environment
env_name <- Sys.getenv("CONDA_DEFAULT_ENV")
if (env_name != "") {
  cat("\nRunning in conda environment:", env_name, "\n")
} else {
  cat("\nNot running in a conda environment\n")
}

cat("\nR version:", R.version.string, "\n")
cat("==================================================\n\n")