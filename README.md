# locusPackRat

A project-based R package for genomic analysis with persistent data storage and Excel output generation.

## Installation

Install directly from GitHub using `devtools`:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install locusPackRat
devtools::install_github("RauLabUNC/locusPackRat")
```

## Core Functions

The package provides three main functions for genomic data analysis:

### 1. `initPackRat()` - Initialize a Project

Create a new analysis project with your gene list or genomic regions:

```r
library(locusPackRat)

# Initialize with a gene list
genes <- data.frame(
  gene_symbol = c("BRCA1", "TP53", "MYC", "EGFR"),
  expression = c(100, 250, 75, 180)
)

initPackRat(
  data = genes,
  mode = "gene",
  species = "human",
  genome = "hg38",
  project_dir = "my_analysis"
)
```

The function automatically:
- Creates a `.locusPackRat` directory structure
- Enriches genes with genomic coordinates from built-in databases (~60,000 mouse genes, ~30,000 human genes)
- Generates orthology mappings between species
- Saves project configuration

### 2. `addRatTable()` - Add Supplementary Data

Link any additional data to your genes or regions:

```r
# Add pathway annotations
pathways <- data.frame(
  gene_symbol = c("BRCA1", "TP53"),
  pathway = c("DNA repair", "Cell cycle"),
  p_value = c(0.001, 0.0001)
)

addRatTable(
  data = pathways,
  table_name = "pathway_analysis",
  link_type = "gene",
  link_by = "gene_symbol",
  project_dir = "my_analysis"
)
```

You can add multiple supplementary tables:
- Expression data
- Phenotype associations
- Clinical information
- QTL statistics
- Any custom data with gene or region identifiers

### 3. `makeGeneSheet()` - Generate Excel Output

Create professional Excel workbooks with all your integrated data:

```r
# Generate Excel with all data
makeGeneSheet(
  format = "excel",
  output_file = "results.xlsx",
  include_supplementary = TRUE,
  project_dir = "my_analysis"
)

# Generate multi-tab Excel with custom filtering
makeGeneSheet(
  format = "excel",
  output_file = "filtered_results.xlsx",
  split_by = "criteria",
  split_criteria = list(
    "High_Expression" = "expression > 100",
    "Low_Expression" = "expression <= 100",
    "Significant" = "p_value < 0.05"
  ),
  highlight_genes = c("BRCA1", "TP53"),
  project_dir = "my_analysis"
)
```

Output features:
- Auto-sized columns with frozen headers
- Gene highlighting (yellow background)
- Multiple tabs based on filtering criteria
- All supplementary data integrated
- Professional formatting ready for publication

## Quick Start Example

```r
library(locusPackRat)
library(data.table)

# 1. Create project with your genes
my_genes <- data.table(
  gene_symbol = c("Myc", "Tp53", "Egfr", "Vegfa", "Il6"),
  fold_change = c(2.5, -1.8, 3.2, 1.5, -2.1),
  p_value = c(0.001, 0.01, 0.0001, 0.05, 0.005)
)

initPackRat(
  data = my_genes,
  mode = "gene",
  species = "mouse",
  genome = "mm39"
)

# 2. Add supplementary data
tissue_exp <- data.table(
  gene_symbol = c("Myc", "Tp53", "Egfr"),
  brain = c(100, 50, 75),
  liver = c(20, 150, 30)
)

addRatTable(
  data = tissue_exp,
  table_name = "tissue_expression",
  link_type = "gene",
  link_by = "gene_symbol"
)

# 3. Generate Excel output
makeGeneSheet(
  format = "excel",
  output_file = "analysis_results.xlsx",
  include_supplementary = TRUE,
  highlight_genes = c("Myc", "Tp53")
)
```

## Features

- **Persistent Storage**: Load data once, use multiple times
- **Built-in Annotations**: Complete genome coordinates for human (hg38) and mouse (mm39)
- **Cross-species Support**: Automatic orthology mapping
- **Flexible Data Integration**: Link any supplementary data by gene symbol or coordinates
- **Professional Excel Output**: Multi-sheet workbooks with formatting and highlighting

## Project Structure

Each project creates a `.locusPackRat` directory:

```
my_analysis/
└── .locusPackRat/
    ├── input/
    │   ├── genes.csv         # Your input genes with coordinates
    │   └── orthology.csv      # Cross-species mappings
    ├── supplementary/
    │   ├── pathway_analysis.csv
    │   └── tissue_expression.csv
    ├── output/
    │   └── results.xlsx
    └── config.json            # Project metadata
```

## Additional Functions

- `listPackRatTables()` - List all supplementary tables in a project
- `generateLocusZoomPlot_v2()` - Create LocusZoom-style visualizations (requires plotgardener)

## Working with Genomic Regions

The package also supports region-based analysis (e.g., QTL intervals, ATAC-seq peaks):

```r
# Initialize with genomic regions
regions <- data.frame(
  chr = c(5, 10, 12),
  start = c(10000000, 50000000, 75000000),
  end = c(15000000, 55000000, 80000000),
  lod_score = c(8.5, 6.2, 7.8)
)

initPackRat(
  data = regions,
  mode = "region",
  species = "mouse",
  genome = "mm39",
  project_dir = "qtl_analysis"
)
```

## Supported Genomes

- Human: hg38 (~30,000 genes)
- Mouse: mm39 (~60,000 genes)

## Dependencies

The package requires:
- R (>= 4.3.0)
- data.table
- jsonlite
- openxlsx
- dplyr
- tidyr
- plotgardener (optional, for visualization)
- RColorBrewer

## Citation

Gural B, Kimball T, Luu A, Rau CD. locusPackRat: A Flexible Framework for Prioritizing Candidate Genes from GWAS and other Gene-Level Studies. *In preparation* (2025).

## Support

- Issues: https://github.com/RauLabUNC/locusPackRat/issues
- Contact: bgural@unc.edu

## License

MIT