# Biomart Gene Coordinate Fetching Scripts

This directory contains scripts for fetching complete gene coordinate information from Ensembl Biomart for various genome assemblies.

## Overview

These scripts query Ensembl Biomart to retrieve gene coordinates (chromosome, start, end, strand) for all genes in specified genome assemblies. The resulting CSV files are used by the genePackRat/locusPackRat package for gene-based analyses.

## Files

### Core Scripts
- `fetch_gene_coordinates.R` - Main R script that queries Biomart
- `submit_all_jobs.sh` - Master script to submit Slurm jobs

### Slurm Job Scripts
- `fetch_mouse_mm39.slurm` - Fetches mouse GRCm39 (mm39) coordinates
- `fetch_human_hg38.slurm` - Fetches human GRCh38 (hg38) coordinates
- `fetch_mouse_mm10.slurm` - Fetches mouse GRCm38 (mm10) coordinates [legacy]
- `fetch_human_hg19.slurm` - Fetches human GRCh37 (hg19) coordinates [legacy]

## Usage

### Submit All Jobs (Recommended)
```bash
# Submit current genome builds only (mm39, hg38)
./submit_all_jobs.sh current

# Submit legacy genome builds only (mm10, hg19)
./submit_all_jobs.sh legacy

# Submit all genome builds
./submit_all_jobs.sh all
```

### Submit Individual Jobs
```bash
# Submit specific genome build
sbatch fetch_mouse_mm39.slurm
sbatch fetch_human_hg38.slurm
```

### Run Locally (Not Recommended for Large Queries)
```bash
# Requires biomaRt and data.table packages
Rscript fetch_gene_coordinates.R mouse mm39
Rscript fetch_gene_coordinates.R human hg38
```

## Output

### File Format
The scripts generate CSV files with the following columns:
- `gene_symbol` - Gene symbol (e.g., BRCA1, Wnt9a)
- `ensembl_id` - Ensembl gene ID (e.g., ENSG00000012048)
- `chr` - Chromosome (1-22/19, X, Y, MT)
- `start` - Gene start position (bp)
- `end` - Gene end position (bp)
- `strand` - Strand orientation (+/-)

### File Locations
Output files are automatically moved to `inst/extdata/`:
- `inst/extdata/mouse_coords_mm39.csv`
- `inst/extdata/human_coords_hg38.csv`
- `inst/extdata/mouse_coords_mm10.csv`
- `inst/extdata/human_coords_hg19.csv`

## Gene Selection Criteria

The scripts filter genes based on:

1. **Chromosomes**: Only standard chromosomes (1-22/19, X, Y, MT)
   - Excludes scaffolds, patches, and alternate assemblies

2. **Gene Biotypes**: Includes biologically relevant types
   - protein_coding
   - lncRNA
   - miRNA, snRNA, snoRNA, rRNA
   - Immunoglobulin genes (IG_*)
   - T-cell receptor genes (TR_*)

3. **Valid Symbols**: Excludes entries without gene symbols

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# View output logs
cat fetch_mm39_coords_*.out

# View error logs (if any)
cat fetch_mm39_coords_*.err
```

## Expected Runtime

- Each job typically completes in 5-15 minutes
- Memory usage: ~2-4 GB per job
- Total genes retrieved:
  - Mouse (mm39/mm10): ~25,000-30,000 genes
  - Human (hg38/hg19): ~20,000-25,000 genes

## Troubleshooting

### Common Issues

1. **Connection timeout**: Ensembl servers may be temporarily unavailable
   - Solution: Wait and resubmit the job

2. **Package installation fails**: BiocManager or biomaRt installation issues
   - Solution: Pre-install packages in an interactive R session

3. **Output directory not found**: inst/extdata doesn't exist
   - Solution: Script automatically creates the directory

### Manual Package Installation
```r
# In an interactive R session
install.packages("BiocManager")
BiocManager::install("biomaRt")
install.packages("data.table")
```

## Ensembl Archive Versions

The scripts use specific Ensembl versions for consistency:
- **mm39**: Current Ensembl release (www.ensembl.org)
- **hg38**: Current Ensembl release (www.ensembl.org)
- **mm10**: Ensembl 109 archive (feb2023.archive.ensembl.org)
- **hg19**: GRCh37 server (grch37.ensembl.org)

## Dependencies

- R >= 4.3.0
- biomaRt (Bioconductor package)
- data.table
- Internet connection to Ensembl servers

## Notes

- The scripts include error handling and progress reporting
- Duplicate genes are automatically removed (keeping first occurrence)
- Results are sorted by chromosome and position
- Email notifications are sent on job completion or failure

## Contact

For issues or questions, contact: bgural@email.unc.edu