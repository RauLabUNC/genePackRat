library(devtools)
library(data.table)
load_all() # Loads your package code

# 1. Initialize a dummy project with one region (Mouse Chr3)
#    (Using the coordinates from your Cisd2 example)
proj_dir <- "test_plot_project"
regions <- data.table(chr = "3", start = 126500000, end = 126600000, region_id = "test_loc")

initPackRat(
  data = regions,
  mode = "region",
  species = "mouse",
  genome = "mm39", 
  project_dir = proj_dir,
  force = TRUE
)

# 2. Create dummy GWAS data
#    plotgardener expects columns: 'chrom', 'pos', 'p' (or 'score')
gwas_data <- data.table(
  chrom = "3",
  pos = seq(126500000, 126600000, by = 2000), # A point every 2kb
  p = 10^-runif(51, 0, 8) # Random p-values (0 to 10^-8)
)

# 3. Save fake data directly to the project
#    (We skip addRatTable() here just to keep the test minimal)
fwrite(gwas_data, file.path(proj_dir, ".locusPackRat/supplementary/test_gwas.csv"))

# 4. Run the plotting function
generateLocusZoomPlot(
  region_id = "test_loc",
  project_dir = proj_dir,
  scan_table = "test_gwas",
  output_file = "test_output.pdf"
)
