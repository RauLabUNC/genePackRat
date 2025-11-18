#!/usr/bin/env Rscript

# Script to prepare test data files for mergeLocusInfo.R functions
# This script creates the necessary test files without running the actual tests

# Set library path for conda environment (miqtl-env)
# This assumes you've run: conda activate miqtl-env
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (conda_prefix != "") {
  conda_lib_path <- file.path(conda_prefix, "lib", "R", "library")
  .libPaths(c(conda_lib_path, .libPaths()))
  cat("Using conda R library path:", conda_lib_path, "\n")
  cat("Full library paths:", paste(.libPaths(), collapse = "\n  "), "\n\n")
} else {
  cat("Warning: CONDA_PREFIX not set.\n")
  cat("Please run: conda activate miqtl-env\n")
  cat("Using default R library paths.\n\n")
}

# Load required libraries
cat("Loading required libraries...\n")
library(data.table)
library(jsonlite)

# Create test data directory
test_data_dir <- "tests/test_data"
if (!dir.exists(test_data_dir)) {
  dir.create(test_data_dir, recursive = TRUE)
  cat("Created test data directory:", test_data_dir, "\n")
}

# ========================================
# PREPARE TEST DATA FILES
# ========================================
cat("\n=== Preparing Test Data Files ===\n")

# 1. Load source data
cat("Loading source data files...\n")
gene_info <- fread("inst/extdata/gene_info.csv")
orthology <- fread("inst/extdata/orthology.csv")
regions <- fread("inst/extdata/all_significant_regions_summary.csv")

cat("  - Loaded", nrow(gene_info), "genes from gene_info.csv\n")
cat("  - Loaded", nrow(orthology), "orthology mappings\n")
cat("  - Loaded", nrow(regions), "regions\n")

# ========================================
# 2. Create Mouse Genome Coordinate Table (mm39)
# ========================================
cat("\nCreating mouse genome coordinate table...\n")

mouse_coords <- gene_info[!is.na(chromosome_name) & chromosome_name != "", .(
  gene_symbol = mouse_gene_symbol,
  ensembl_id = mouse_ensembl_id,
  chr = chromosome_name,
  start = start_position,
  end = end_position,
  strand = "+"  # Placeholder - we don't have strand info
)]

# Clean up: remove NA gene symbols and duplicates
mouse_coords <- mouse_coords[!is.na(gene_symbol) & gene_symbol != ""]
mouse_coords <- unique(mouse_coords, by = "gene_symbol")

# Save mouse coordinates
mouse_coords_file <- file.path(test_data_dir, "mouse_coords_mm39.csv")
fwrite(mouse_coords, mouse_coords_file)
cat("  Saved", nrow(mouse_coords), "mouse genes to", mouse_coords_file, "\n")

# ========================================
# 3. Create Human Coordinate Table (hg38)
# ========================================
cat("\nCreating human genome coordinate table...\n")

# Get unique human genes from orthology
human_genes <- unique(orthology[!is.na(human_gene_symbol) & human_gene_symbol != "", .(
  gene_symbol = human_gene_symbol,
  ensembl_id = human_ensembl_id
)])

# For testing, create realistic dummy coordinates for human genes
# We'll base them on the mouse coordinates where we have orthologs
human_coords <- merge(
  human_genes,
  merge(orthology, mouse_coords, by.x = "mouse_ensembl_id", by.y = "ensembl_id"),
  by.x = c("gene_symbol", "ensembl_id"),
  by.y = c("human_gene_symbol", "human_ensembl_id")
)[, .(
  gene_symbol = gene_symbol,
  ensembl_id = ensembl_id,
  chr = chr,  # Use same chromosome for simplicity in testing
  start = start + sample(-1000:1000, .N, replace = TRUE),  # Slightly offset positions
  end = end + sample(-1000:1000, .N, replace = TRUE),
  strand = "+"
)]

# Save human coordinates
human_coords_file <- file.path(test_data_dir, "human_coords_hg38.csv")
fwrite(human_coords, human_coords_file)
cat("  Saved", nrow(human_coords), "human genes to", human_coords_file, "\n")

# ========================================
# 4. Create Orthology Table
# ========================================
cat("\nCreating orthology table...\n")

ortho_table <- merge(
  orthology[, .(mouse_ensembl_id, human_ensembl_id, human_gene_symbol)],
  gene_info[, .(mouse_ensembl_id, mouse_gene_symbol)],
  by = "mouse_ensembl_id",
  all.x = TRUE
)

# Save orthology table
ortho_file <- file.path(test_data_dir, "orthology_mapping.csv")
fwrite(ortho_table, ortho_file)
cat("  Saved", nrow(ortho_table), "orthology mappings to", ortho_file, "\n")

# ========================================
# 5. Create Test Gene List
# ========================================
cat("\nCreating test gene list...\n")

# Sample 100 random mouse genes for testing
set.seed(42)
test_genes <- mouse_coords[sample(.N, min(100, .N)), .(
  gene_symbol = gene_symbol,
  custom_value = rnorm(.N),
  category = sample(c("A", "B", "C"), .N, replace = TRUE)
)]

# Save test gene list
test_genes_file <- file.path(test_data_dir, "test_gene_list.csv")
fwrite(test_genes, test_genes_file)
cat("  Saved", nrow(test_genes), "test genes to", test_genes_file, "\n")

# ========================================
# 6. Create Test Region Data (Chr 5)
# ========================================
cat("\nCreating test region data...\n")

# Filter to chromosome 5 and convert column names
test_regions <- regions[chr == "5", .(
  chr = chr,
  start = upper_pos_lod_drop * 1e6,  # Convert Mb to bp
  end = lower_pos_lod_drop * 1e6,
  peak_pos = peak_pos,
  max_lod = max_lod,
  trait = trait,
  drug = drug
)]

# Ensure start < end (swap if needed)
swap_idx <- test_regions$start > test_regions$end
if (any(swap_idx)) {
  temp_start <- test_regions[swap_idx, start]
  test_regions[swap_idx, start := end]
  test_regions[swap_idx, end := temp_start]
}

# Add region IDs
test_regions[, region_id := paste0("region_chr5_", .I)]

# Save test regions
test_regions_file <- file.path(test_data_dir, "test_regions_chr5.csv")
fwrite(test_regions, test_regions_file)
cat("  Saved", nrow(test_regions), "chromosome 5 regions to", test_regions_file, "\n")

# ========================================
# 7. Create Test Expression Data
# ========================================
cat("\nCreating test expression data...\n")

# Create expression data for subset of test genes
expression_data <- data.table(
  gene_symbol = test_genes$gene_symbol[1:min(50, nrow(test_genes))],
  expression_ctrl = runif(min(50, nrow(test_genes)), 0, 100),
  expression_treated = runif(min(50, nrow(test_genes)), 0, 100),
  log2fc = rnorm(min(50, nrow(test_genes))),
  pvalue = runif(min(50, nrow(test_genes)), 0, 1)
)

# Save expression data
expression_file <- file.path(test_data_dir, "test_expression_data.csv")
fwrite(expression_data, expression_file)
cat("  Saved expression data for", nrow(expression_data), "genes to", expression_file, "\n")

# ========================================
# 8. Create Test Phenotype Data
# ========================================
cat("\nCreating test phenotype data...\n")

# Create phenotype associations for some test genes
phenotype_data <- data.table(
  gene_symbol = rep(test_genes$gene_symbol[1:20], each = 2),
  phenotype = rep(c("cardiac dysfunction", "abnormal heart rate",
                    "increased body weight", "decreased lifespan"), 10),
  p_value = runif(40, 0, 0.1),
  source = "MGI"
)

# Save phenotype data
phenotype_file <- file.path(test_data_dir, "test_phenotype_data.csv")
fwrite(phenotype_data, phenotype_file)
cat("  Saved phenotype data to", phenotype_file, "\n")

# ========================================
# SUMMARY
# ========================================
cat("\n=== Test Data Preparation Complete ===\n")
cat("Test data directory:", test_data_dir, "\n")
cat("Files created:\n")
test_files <- list.files(test_data_dir, full.names = FALSE)
for (f in test_files) {
  file_size <- file.info(file.path(test_data_dir, f))$size
  cat(sprintf("  - %s (%.1f KB)\n", f, file_size/1024))
}

cat("\nNext steps:\n")
cat("1. Review the test data files in", test_data_dir, "\n")
cat("2. Run the actual tests with: Rscript tests/run_tests.R\n")
cat("3. Test data can be used to validate the three main functions:\n")
cat("   - initPackRat() with test_gene_list.csv or test_regions_chr5.csv\n")
cat("   - addRatTable() with test_expression_data.csv or test_phenotype_data.csv\n")
cat("   - generateGeneSheet() after initialization\n")