#!/bin/bash
#SBATCH --job-name=fetch_hg19_coords
#SBATCH --output=scripts/biomart/logs/fetch_hg19_coords_%j.out
#SBATCH --error=scripts/biomart/logs/fetch_hg19_coords_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8GB
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Print job info
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo "Working Directory: $(pwd)"
echo "=========================================="
echo ""

# Load R module
module load r/4.4.0
echo "R version:"
R --version | head -n 1
echo ""

# Create output directory if it doesn't exist
OUTPUT_DIR="${PWD}/inst/extdata"
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# Install/check biomaRt package
echo "Checking/Installing biomaRt package..."
Rscript -e '
if (!require("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
  BiocManager::install("biomaRt", ask = FALSE, quiet = TRUE)
}
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org", quiet = TRUE)
}
cat("biomaRt version:", as.character(packageVersion("biomaRt")), "\n")
cat("data.table version:", as.character(packageVersion("data.table")), "\n")
'
echo ""

# Change to script directory
cd scripts/biomart

# Run the R script to fetch human hg19 coordinates
echo "=========================================="
echo "Fetching human hg19 gene coordinates..."
echo "=========================================="
Rscript fetch_gene_coordinates.R human hg19

# Check if the script succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "Successfully generated human_coords_hg19.csv"

    # Move the file to inst/extdata
    if [ -f "human_coords_hg19.csv" ]; then
        mv human_coords_hg19.csv "$OUTPUT_DIR/"
        echo "Moved file to: $OUTPUT_DIR/human_coords_hg19.csv"

        # Print file info
        echo ""
        echo "File statistics:"
        wc -l "$OUTPUT_DIR/human_coords_hg19.csv"
        echo ""
        echo "First 5 lines:"
        head -n 5 "$OUTPUT_DIR/human_coords_hg19.csv"
    else
        echo "ERROR: Output file not found!"
        exit 1
    fi
else
    echo "ERROR: R script failed!"
    exit 1
fi

echo ""
echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="