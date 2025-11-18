# Test installation script
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos="https://cloud.r-project.org")
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos="https://cloud.r-project.org")
}

# Install locusPackRat locally
remotes::install_local("/proj/raulab/users/brian/packrat", force = TRUE, upgrade = "never")

# Test loading the package
library(locusPackRat)

# Check that functions are available
cat("Testing function availability:\n")
cat("  initPackRat:", exists("initPackRat"), "\n")
cat("  addRatTable:", exists("addRatTable"), "\n")
cat("  makeGeneSheet:", exists("makeGeneSheet"), "\n")
cat("  listPackRatTables:", exists("listPackRatTables"), "\n")
cat("  filterGenes:", exists("filterGenes"), "\n")
cat("  makeFilter:", exists("makeFilter"), "\n")

# Test with embedded coordinate files
cat("\nTesting embedded coordinate files:\n")
mouse_coords <- system.file("extdata", "mouse_coords_mm39.csv", package = "locusPackRat")
human_coords <- system.file("extdata", "human_coords_hg38.csv", package = "locusPackRat")
cat("  Mouse coords found:", file.exists(mouse_coords), "\n")
cat("  Human coords found:", file.exists(human_coords), "\n")

cat("\nPackage installation test complete!\n")