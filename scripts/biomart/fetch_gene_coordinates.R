#!/usr/bin/env Rscript

#' Fetch Gene Coordinates from Ensembl Biomart
#'
#' This script queries Ensembl Biomart to retrieve complete gene coordinates
#' for specified genome assemblies. It generates CSV files with gene symbols,
#' Ensembl IDs, and genomic coordinates.
#'
#' Usage:
#'   Rscript fetch_gene_coordinates.R <species> <genome>
#'
#' Arguments:
#'   species: "mouse" or "human"
#'   genome:  "mm39", "mm10", "hg38", or "hg19"
#'
#' Output:
#'   Creates a CSV file: <species>_coords_<genome>.csv

# Load required libraries
suppressPackageStartupMessages({
  library(biomaRt)
  library(data.table)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Usage: Rscript fetch_gene_coordinates.R <species> <genome>\n")
  cat("  species: mouse or human\n")
  cat("  genome:  mm39, mm10, hg38, or hg19\n")
  quit(status = 1)
}

species <- tolower(args[1])
genome <- tolower(args[2])

# Validate arguments
valid_combinations <- list(
  mouse = c("mm39", "mm10"),
  human = c("hg38", "hg19")
)

if (!species %in% names(valid_combinations)) {
  stop(paste("Invalid species:", species, ". Must be 'mouse' or 'human'"))
}

if (!genome %in% valid_combinations[[species]]) {
  stop(paste("Invalid genome for", species, ":", genome,
             ". Must be one of:", paste(valid_combinations[[species]], collapse = ", ")))
}

# Function to get Ensembl dataset and version
get_ensembl_info <- function(species, genome) {
  info <- list()

  if (species == "mouse") {
    info$dataset <- "mmusculus_gene_ensembl"
    info$species_name <- "Mus musculus"

    if (genome == "mm39") {
      # mm39 corresponds to GRCm39 (Ensembl release 104+)
      info$host <- "https://www.ensembl.org"
      info$version <- "GRCm39"
    } else if (genome == "mm10") {
      # mm10 corresponds to GRCm38 (use archive for stable version)
      info$host <- "https://feb2023.archive.ensembl.org"  # Ensembl 109
      info$version <- "GRCm38"
    }
  } else if (species == "human") {
    info$dataset <- "hsapiens_gene_ensembl"
    info$species_name <- "Homo sapiens"

    if (genome == "hg38") {
      # hg38 corresponds to GRCh38
      info$host <- "https://www.ensembl.org"
      info$version <- "GRCh38"
    } else if (genome == "hg19") {
      # hg19 corresponds to GRCh37 (use dedicated server)
      info$host <- "https://grch37.ensembl.org"
      info$version <- "GRCh37"
    }
  }

  return(info)
}

# Get Ensembl configuration
ensembl_info <- get_ensembl_info(species, genome)

cat("\n========================================\n")
cat("Fetching gene coordinates from Ensembl\n")
cat("========================================\n")
cat("Species:", ensembl_info$species_name, "\n")
cat("Genome assembly:", genome, "(", ensembl_info$version, ")\n")
cat("Dataset:", ensembl_info$dataset, "\n")
cat("Host:", ensembl_info$host, "\n")
cat("========================================\n\n")

# Clear BiocFileCache to avoid cache corruption issues
cat("Clearing BiocFileCache to avoid corruption issues...\n")
tryCatch({
  if (requireNamespace("BiocFileCache", quietly = TRUE)) {
    cache_dir <- tools::R_user_dir("biomaRt", "cache")
    if (dir.exists(cache_dir)) {
      unlink(cache_dir, recursive = TRUE)
      cat("Cache cleared successfully\n")
    }
  }
}, error = function(e) {
  cat("Note: Could not clear cache, continuing anyway...\n")
})

# Connect to Ensembl
cat("Connecting to Ensembl Biomart...\n")
tryCatch({
  mart <- useMart(
    biomart = "ensembl",
    dataset = ensembl_info$dataset,
    host = ensembl_info$host
  )
}, error = function(e) {
  cat("Error connecting to Biomart:\n")
  cat(e$message, "\n")
  quit(status = 1)
})

cat("Connection established!\n\n")

# Define attributes to retrieve
if (species == "mouse") {
  attributes <- c(
    "mgi_symbol",           # Mouse gene symbol
    "ensembl_gene_id",      # Ensembl gene ID
    "chromosome_name",      # Chromosome
    "start_position",       # Gene start
    "end_position",         # Gene end
    "strand",               # Strand (1 or -1)
    "gene_biotype"          # Gene biotype for filtering
  )
} else {
  attributes <- c(
    "hgnc_symbol",          # Human gene symbol
    "ensembl_gene_id",      # Ensembl gene ID
    "chromosome_name",      # Chromosome
    "start_position",       # Gene start
    "end_position",         # Gene end
    "strand",               # Strand (1 or -1)
    "gene_biotype"          # Gene biotype for filtering
  )
}

# Query Biomart
cat("Querying Biomart for gene coordinates...\n")
cat("This may take several minutes...\n\n")

start_time <- Sys.time()
tryCatch({
  genes <- getBM(
    attributes = attributes,
    mart = mart
  )
}, error = function(e) {
  cat("Error querying Biomart:\n")
  cat(e$message, "\n")
  quit(status = 1)
})
end_time <- Sys.time()

cat("Query completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
cat("Retrieved", nrow(genes), "gene entries\n\n")

# Convert to data.table for efficient processing
setDT(genes)

# Standardize column names
if (species == "mouse") {
  setnames(genes, "mgi_symbol", "gene_symbol")
} else {
  setnames(genes, "hgnc_symbol", "gene_symbol")
}

# Filter for standard chromosomes only (exclude scaffolds, patches, etc.)
cat("Filtering for standard chromosomes...\n")
standard_chr <- c(as.character(1:50), "X", "Y", "MT")  # Cover human (1-22) and mouse (1-19)
genes <- genes[chromosome_name %in% standard_chr]
cat("Genes on standard chromosomes:", nrow(genes), "\n\n")

# Filter for protein-coding genes and other important biotypes
cat("Filtering by gene biotype...\n")
important_biotypes <- c(
  "protein_coding",
  "lncRNA",
  "miRNA",
  "snRNA",
  "snoRNA",
  "rRNA",
  "IG_C_gene",
  "IG_D_gene",
  "IG_J_gene",
  "IG_V_gene",
  "TR_C_gene",
  "TR_D_gene",
  "TR_J_gene",
  "TR_V_gene"
)

genes_filtered <- genes[gene_biotype %in% important_biotypes]
cat("Genes after biotype filtering:", nrow(genes_filtered), "\n")
cat("Biotypes included:", paste(unique(genes_filtered$gene_biotype), collapse = ", "), "\n\n")

# Remove genes without symbols
cat("Removing entries without gene symbols...\n")
genes_filtered <- genes_filtered[!is.na(gene_symbol) & gene_symbol != ""]
cat("Genes with valid symbols:", nrow(genes_filtered), "\n\n")

# Convert strand from numeric to +/- format
genes_filtered[, strand := ifelse(strand == 1, "+", "-")]

# Select and order final columns
final_genes <- genes_filtered[, .(
  gene_symbol,
  ensembl_id = ensembl_gene_id,
  chr = chromosome_name,
  start = start_position,
  end = end_position,
  strand
)]

# Remove duplicates (keeping first occurrence)
cat("Removing duplicate entries...\n")
initial_count <- nrow(final_genes)
final_genes <- unique(final_genes, by = "ensembl_id")
duplicates_removed <- initial_count - nrow(final_genes)
cat("Duplicates removed:", duplicates_removed, "\n")
cat("Final gene count:", nrow(final_genes), "\n\n")

# Sort by chromosome and position
# Custom chromosome ordering
chr_order <- c(as.character(1:50), "X", "Y", "MT")
final_genes[, chr_factor := factor(chr, levels = chr_order)]
setorder(final_genes, chr_factor, start)
final_genes[, chr_factor := NULL]

# Generate output filename
output_file <- paste0(species, "_coords_", genome, ".csv")
output_path <- file.path(getwd(), output_file)

# Save to CSV
cat("Saving results to:", output_path, "\n")
fwrite(final_genes, output_path, row.names = FALSE)

# Print summary statistics
cat("\n========================================\n")
cat("Summary Statistics\n")
cat("========================================\n")
cat("Total genes:", nrow(final_genes), "\n")
cat("Unique chromosomes:", length(unique(final_genes$chr)), "\n")
cat("Chromosomes:", paste(sort(unique(final_genes$chr)), collapse = ", "), "\n")
cat("Genes on plus strand:", sum(final_genes$strand == "+"), "\n")
cat("Genes on minus strand:", sum(final_genes$strand == "-"), "\n")

# Sample of the data
cat("\nFirst 5 entries:\n")
print(head(final_genes, 5))

cat("\n========================================\n")
cat("Gene coordinate file created successfully!\n")
cat("Output:", output_path, "\n")
cat("========================================\n\n")