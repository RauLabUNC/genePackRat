# Test script for makeGeneSheet features
# Tests: overlap_mode, overlap_type column, exclude_tables, removeRatTable

devtools::load_all()
library(data.table)

# Clean up any previous test directories
if (dir.exists("test_region_project")) unlink("test_region_project", recursive = TRUE)
if (dir.exists("test_complete_overlap")) unlink("test_complete_overlap", recursive = TRUE)

# ============================================================
# 1. Create a region-mode project with default overlap_mode="any"
# ============================================================
message("\n=== Creating test region project (overlap_mode='any') ===")

# Define two test regions on chr1 (gene-dense areas)
# Region 1: chr1:3.2M-4M
#   - Xkr4 (3276124-3741721) is COMPLETE (entirely within region)
#   - Gm71132 (3943429-4018627) is PARTIAL (extends past 4M)
#   - Gm37381 (3975962-4056438) is PARTIAL (extends past 4M)
# Region 2: chr1:4.8M-5.7M - contains several complete genes
test_regions <- data.table(
  region_id = c("region_1", "region_2"),
  chr = c("1", "1"),
  start = c(3200000, 4800000),
  end = c(4000000, 5700000)
)

initPackRat(
  data = test_regions,
  mode = "region",
  species = "mouse",
  genome = "mm39",
  project_dir = "test_region_project",
  overlap_mode = "any",  # default - keep partial and complete overlaps
  keep_pseudo = TRUE,    # include Gm* genes so we can see partial overlaps
  force = TRUE
)

# ============================================================
# 2. Verify overlap_types in regions.csv
# ============================================================
message("\n=== Verifying overlap_types in stored regions ===")

regions <- fread("test_region_project/.locusPackRat/input/regions.csv")
message("Region 1:")
message("  Genes: ", regions$genes[1])
message("  Overlap types: ", regions$overlap_types[1])
message("Region 2:")
message("  Genes: ", regions$genes[2])
message("  Overlap types: ", regions$overlap_types[2])

# ============================================================
# 3. Add supplementary tables for testing
# ============================================================
message("\n=== Adding supplementary tables ===")

# 3a. SNP data with genomic coordinates (tests interval overlap with genes)
snp_data <- data.table(
  chr = c("1", "1", "1", "1", "1", "1"),
  start = c(3300000, 3550000, 3850000, 4880000, 5050000, 5660000),
  end = c(3300000, 3550000, 3850000, 4880000, 5050000, 5660000),
  snp_id = c("rs001", "rs002", "rs003", "rs004", "rs005", "rs006"),
  effect = c("missense", "synonymous", "intergenic", "intronic", "missense", "stop")
)

addRatTable(
  data = snp_data,
  table_name = "test_snps",
  abbreviation = "snp",
  link_type = "point",
  project_dir = "test_region_project"
)

# 3b. Gene-level annotation
genes_in_regions <- unlist(strsplit(regions$genes, ", "))
genes_in_regions <- genes_in_regions[!is.na(genes_in_regions)]
gene_annotations <- data.table(
  gene_symbol = genes_in_regions,
  score = runif(length(genes_in_regions), 0, 1),
  category = sample(c("A", "B", "C"), length(genes_in_regions), replace = TRUE)
)

addRatTable(
  data = gene_annotations,
  table_name = "gene_scores",
  abbreviation = "gs",
  link_type = "gene",
  project_dir = "test_region_project"
)

# 3c. Founder data (for exclude_tables test)
founder_data <- data.table(
  chr = c("1", "1"),
  start = c(3200000, 4800000),
  end = c(4000000, 5700000),
  lod = c(5.2, 3.8),
  strain_pattern = c("ABCDEFGH", "HGFEDCBA")
)

addRatTable(
  data = founder_data,
  table_name = "founder_info",
  abbreviation = "fi",
  link_type = "region",
  project_dir = "test_region_project"
)

# ============================================================
# 4. Test makeGeneSheet - verify overlap_type column
# ============================================================
message("\n=== Testing makeGeneSheet output ===")

result <- makeGeneSheet(
  output_file = "test_default.xlsx",
  project_dir = "test_region_project"
)

message("Output columns: ", paste(names(result), collapse = ", "))
message("Number of genes: ", nrow(result))

# Check overlap_type column
message("\n--- overlap_type column ---")
if ("overlap_type" %in% names(result)) {
  message("SUCCESS: overlap_type column present")
  print(result[, .(gene_symbol, overlap_type)])
} else {
  message("ERROR: overlap_type column missing!")
}

# Check SNP assignments (should be gene-specific, not region-wide)
message("\n--- SNP assignments (should be gene-specific) ---")
snp_col <- intersect(c("snp_snp_id", "snp_id"), names(result))[1]
if (!is.na(snp_col)) {
  print(result[, .(gene_symbol, snps = get(snp_col))])
}

# ============================================================
# 5. Test exclude_tables parameter
# ============================================================
message("\n=== Testing exclude_tables parameter ===")

result_excluded <- makeGeneSheet(
  output_file = "test_excluded.xlsx",
  exclude_tables = c("founder_info"),
  project_dir = "test_region_project"
)

founder_cols <- grep("founder|fi_|lod|strain", names(result_excluded), value = TRUE)
if (length(founder_cols) == 0) {
  message("SUCCESS: founder_info columns excluded")
} else {
  message("ERROR: founder columns still present: ", paste(founder_cols, collapse = ", "))
}

# ============================================================
# 6. Test removeRatTable function
# ============================================================
message("\n=== Testing removeRatTable ===")

message("Before removal:")
listPackRatTables("test_region_project")

removeRatTable("gene_scores", project_dir = "test_region_project")

message("\nAfter removal:")
listPackRatTables("test_region_project")

# ============================================================
# 7. Test overlap_mode = "complete"
# ============================================================
message("\n=== Testing overlap_mode='complete' ===")

initPackRat(
  data = test_regions,
  mode = "region",
  species = "mouse",
  genome = "mm39",
  project_dir = "test_complete_overlap",
  overlap_mode = "complete",  # only keep genes completely within regions
  keep_pseudo = TRUE,         # include Gm* genes for comparison
  force = TRUE
)

regions_complete <- fread("test_complete_overlap/.locusPackRat/input/regions.csv")

message("With overlap_mode='complete':")
message("  Region 1 genes: ", regions_complete$genes[1])
message("  Region 1 types: ", regions_complete$overlap_types[1])
message("  Region 2 genes: ", regions_complete$genes[2])
message("  Region 2 types: ", regions_complete$overlap_types[2])

# Compare gene counts
count_genes <- function(genes_col) {
  sum(sapply(strsplit(genes_col, ", "), function(x) {
    if (length(x) == 1 && is.na(x)) 0 else length(x)
  }))
}

any_count <- count_genes(regions$genes)
complete_count <- count_genes(regions_complete$genes)

message("\nGene count comparison:")
message("  overlap_mode='any': ", any_count, " genes")
message("  overlap_mode='complete': ", complete_count, " genes")

if (complete_count < any_count) {
  message("SUCCESS: 'complete' mode has fewer genes (partial overlaps filtered)")
} else if (complete_count == any_count) {
  message("NOTE: All genes have complete overlaps in these regions")
}

# ============================================================
# 8. Summary
# ============================================================
message("\n=== Test Summary ===")
message("1. overlap_mode parameter: TESTED")
message("2. overlap_type column in output: ",
        ifelse("overlap_type" %in% names(result), "PASS", "FAIL"))
message("3. exclude_tables parameter: ",
        ifelse(length(founder_cols) == 0, "PASS", "FAIL"))
message("4. removeRatTable function: TESTED")
message("5. Gene-specific SNP overlap: TESTED")

message("\nOutput files in: test_region_project/.locusPackRat/output/")
