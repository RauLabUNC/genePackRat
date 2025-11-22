library(devtools)
library(data.table)
load_all() 

# ========================================================
# 1. Load the Real Data
# ========================================================
scan_file <- "EF.21_Iso.rds"
thresh_file <- "EF.21_Iso_threshold.rds"

if (!file.exists(scan_file)) stop("Could not find ", scan_file)

message("Loading scan data from ", scan_file, "...")
real_scan <- readRDS(scan_file)

# ========================================================
# 2. Extract Data from Complex List Structure
# ========================================================
message("Extracting data vectors...")

# Create base table from the list elements
# We use standard names 'chrom' and 'pos' which addRatTable will recognize
gwas_dt <- data.table(
  chrom  = as.character(real_scan$chr),           # Chromosomes
  pos    = as.numeric(real_scan$pos$Mb) * 1e6,    # Convert Mb to Basepairs
  p      = as.numeric(real_scan$p.value),         # P-values
  lod    = as.numeric(real_scan$LOD),             # LOD scores
  marker = as.character(real_scan$loci)           # Marker IDs
)

# OPTIONAL: Convert LOD to pseudo-p-value for plotting
# (plotgardener plots -log10(p), so using 10^-LOD makes the Y-axis equal the LOD score)
if (TRUE) {
  message("Converting LOD to pseudo-p-values for plotting...")
  gwas_dt[, p := 10^-lod]
}

# ========================================================
# 3. Define 2MB Region (Chr 3 Peak) & Subset
# ========================================================
target_chrom <- "3"
center_pos   <- 138860000
start_pos    <- center_pos - 1000000 # 1MB left
end_pos      <- center_pos + 1000000 # 1MB right
region_id    <- "EF21_Iso_Peak_Zoom"

message(sprintf("Subsetting for Chr %s: %.2fMb - %.2fMb...", 
                target_chrom, start_pos/1e6, end_pos/1e6))

# Subset the GWAS data
# Note: We only need to filter rows. We don't need to change columns anymore.
subset_gwas <- gwas_dt[chrom == target_chrom & pos >= start_pos & pos <= end_pos]

if (nrow(subset_gwas) == 0) stop("No SNPs found in the requested region!")

message("Found ", nrow(subset_gwas), " markers in the target window.")

# ========================================================
# 4. Initialize Project & Plot
# ========================================================
proj_dir <- "EF21_Iso_Analysis"

# Initialize the project
initPackRat(
  data = data.table(chr = target_chrom, start = start_pos, end = end_pos, region_id = region_id),
  mode = "region",
  species = "mouse",
  genome = "mm39",
  project_dir = proj_dir,
  force = TRUE
)

# addRatTable will automatically find 'chrom' and 'pos', create intervals, and link them.
addRatTable(
  data = gwas_dt,
  table_name = "full_iso_scan",
  link_type = "point", 
  project_dir = proj_dir
)

# Generate Plot
pdf_name <- "EF21_Iso_LocusZoom.pdf"
generateLocusZoomPlot(
  region_id = region_id,
  project_dir = proj_dir,
  scan_table = "full_iso_scan",
  output_file = pdf_name
)
