# locusPackRat

An R package for managing genomic analysis projects. Stores genes or regions with supplementary data and generates formatted Excel output.

## Installation

```r
devtools::install_github("RauLabUNC/locusPackRat")
```

## Quick Start

```r
library(locusPackRat)

# 1. Initialize a project with your genes
genes <- data.frame(
  gene_symbol = c("Myc", "Tp53", "Egfr", "Vegfa"),
  log2FC = c(2.5, -1.8, 3.2, 1.5),
  padj = c(0.001, 0.01, 0.0001, 0.05)
)

initPackRat(
  data = genes,
  mode = "gene",
  species = "mouse",
  genome = "mm39",
  project_dir = "my_analysis"
)

# 2. Query external databases (optional)
queryMouseMine(project_dir = "my_analysis")
queryOpenTargets(project_dir = "my_analysis")

# 3. Generate Excel output
makeGeneSheet(
  format = "excel",
  output_file = "results.xlsx",
  project_dir = "my_analysis"
)
```

## Core Functions

### `initPackRat()`

Creates a project directory with your genes or genomic regions. Automatically adds genomic coordinates and orthology mappings from built-in reference files.

```r
initPackRat(
  data = my_genes,
  mode = "gene",
  species = "mouse",
  genome = "mm39",
  project_dir = "."
)
```

**Modes:** `"gene"` or `"region"`
**Species:** `"human"` or `"mouse"`
**Genomes:** `"hg38"`, `"hg19"`, `"mm39"`, `"mm10"`

### `addRatTable()`

Links supplementary data to your project by gene symbol or coordinates.

```r
addRatTable(
  data = my_annotations,
  table_name = "pathway_data",
  link_type = "gene",
  link_by = "gene_symbol",
  abbreviation = "pw",
  project_dir = "."
)
```

The `abbreviation` parameter sets a short prefix (e.g., "pw") used to distinguish columns when multiple tables share the same column names in Excel output.

### `makeGeneSheet()`

Generates CSV or Excel output with all project data merged together.

```r
makeGeneSheet(
  format = "excel",
  output_file = "results.xlsx",
  include_supplementary = TRUE,
  prefix_mode = "collision",
  project_dir = "."
)
```

**prefix_mode options:**
- `"collision"` (default): Only prefixes columns that would overwrite existing ones
- `"abbreviated"`: Prefixes all columns from tables that have an abbreviation set
- `"always"`: Prefixes all supplementary columns using abbreviation or table name

## Query Functions

These functions fetch data from external databases and save results as supplementary tables.

### `queryMouseMine()`

Queries MouseMine for gene-phenotype associations (mouse projects only). Uses the REST API directly—no InterMineR dependency.

```r
queryMouseMine(project_dir = "my_analysis")
# Creates: mouse_phenotypes table (abbreviation: "mm")
```

### `queryOpenTargets()`
Queries Open Targets for disease associations, genetic constraints, and drug tractability. Works for both human and mouse projects (mouse genes are mapped to human orthologs).

```r
queryOpenTargets(project_dir = "my_analysis")
# Creates: ot_diseases (abbrev: "otd"), ot_constraints ("otc"), ot_tractability ("ott")
```

## Visualization Function

###`generateLocusZoomPlot()`
Generates a LocusZoom-style plot for visualization of a locus.  Requires the plotGardener package.  

```r 
generateLocusZoomPlot(
     region_id="region_1",
     project_dir="my_analysis",
     scan_table="LoD_Values",
     signal_table="CC_Founder_Data",
     width=10,
     height=6,
     threshold=4,
     layout_ratios = c(manhattan = 0.35, signal = 0.40, genes = 0.25))

```


## Other Functions

- `listPackRatTables()` - List supplementary tables in a project
- `removeRatTable()` - Delete a supplementary table from a project


## Project Structure

```
my_analysis/
└── .locusPackRat/
    ├── input/
    │   ├── genes.csv
    │   └── orthology.csv
    ├── supplementary/
    │   ├── mouse_phenotypes.csv
    │   └── ot_diseases.csv
    ├── output/
    │   └── results.xlsx
    └── config.json
```

## Supported Genomes

| Species | Genome | Genes |
|---------|--------|-------|
| Human | hg38 | ~30,000 |
| Human | hg19 | ~24,000 |
| Mouse | mm39 | ~57,000 |
| Mouse | mm10 | ~38,000 |

## Dependencies

**Required:** data.table, jsonlite, openxlsx, dplyr, httr, Bioconductor 3.22
**Optional:** plotgardener (for visualization)

## Citation

Gural B, Kimball T, Luu A, Rau CD. locusPackRat: A Flexible Framework for Prioritizing Candidate Genes from GWAS and other Gene-Level Studies. *Under Review* (2025).

## License

MIT
