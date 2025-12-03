# Tests for generateUpsetPlot, generateLocusZoomPlot, and checkCriteria

test_that("generateUpsetPlot creates output with minimal setup", {
  test_dir <- file.path(tempdir(), "test_upset_minimal")
  on.exit(unlink(test_dir, recursive = TRUE))

  # Create simple gene-mode project
  test_genes <- c("Myc", "Tp53", "Egfr", "Vegfa", "Gapdh", "Brca1", "Pten", "Kras")
  test_data <- data.frame(
    gene_symbol = test_genes,
    log2FC = c(1.5, -2.1, 0.8, 1.2, -0.3, 2.0, -1.5, 0.5)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Add two supplementary tables
  set1_data <- data.table::data.table(
    gene_symbol = c("Myc", "Tp53", "Egfr", "Kras"),
    score = c(0.01, 0.02, 0.5, 0.03)
  )
  data.table::fwrite(
    set1_data,
    file.path(test_dir, ".locusPackRat", "supplementary", "set1.csv")
  )

  set2_data <- data.table::data.table(
    gene_symbol = c("Tp53", "Vegfa", "Brca1", "Kras"),
    category = c("A", "B", "A", "C")
  )
  data.table::fwrite(
    set2_data,
    file.path(test_dir, ".locusPackRat", "supplementary", "set2.csv")
  )

  criteria <- list(
    "Low Score" = list(table = "set1", column = "score", condition = "< 0.05"),
    "In Set 2" = list(table = "set2", column = "gene_symbol", condition = "exists")
  )

  result <- generateUpsetPlot(
    criteria = criteria,
    project_dir = test_dir,
    output_file = "minimal_upset.pdf"
  )

  expect_true(file.exists(
    file.path(test_dir, ".locusPackRat", "output", "minimal_upset.pdf")
  ))
})

test_that("generateLocusZoomPlot creates output with minimal setup", {
  test_dir <- file.path(tempdir(), "test_zoom_minimal")
  on.exit(unlink(test_dir, recursive = TRUE))

  # Create region-mode project
  region_data <- data.frame(
    chr = "1",
    start = 1000000,
    end = 2000000,
    region_id = "test_region"
  )

  initPackRat(
    data = region_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Add scan data
  positions <- seq(1000000, 2000000, by = 10000)
  scan_data <- data.table::data.table(
    chr = "1",
    pos = positions,
    start = positions,
    end = positions,
    lod = abs(rnorm(length(positions), 2, 1))
  )
  # Add a peak
  scan_data[50, lod := 5.5]

  addRatTable(
    data = scan_data,
    table_name = "test_scan",
    abbreviation = "ts",
    link_type = "region",
    link_by = "chr,start,end",
    project_dir = test_dir
  )

  result <- generateLocusZoomPlot(
    region_id = "test_region",
    project_dir = test_dir,
    scan_table = "test_scan",
    threshold = 4,
    output_file = "minimal_zoom.pdf"
  )

  expect_true(file.exists(
    file.path(test_dir, ".locusPackRat", "output", "minimal_zoom.pdf")
  ))
})

test_that("checkCriteria validates correct criteria", {
  test_dir <- file.path(tempdir(), "test_check_valid")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53", "Egfr"))

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  supp_data <- data.table::data.table(
    gene_symbol = c("Myc", "Tp53"),
    value = c(0.01, 0.5)
  )
  data.table::fwrite(
    supp_data,
    file.path(test_dir, ".locusPackRat", "supplementary", "test_table.csv")
  )

  criteria <- list(
    "Low Value" = list(table = "test_table", column = "value", condition = "< 0.1")
  )

  expect_true(checkCriteria(criteria, project_dir = test_dir, verbose = FALSE))
})

test_that("checkCriteria errors on missing table", {
  test_dir <- file.path(tempdir(), "test_check_missing")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  criteria <- list(
    "Bad" = list(table = "nonexistent_table", column = "x", condition = "< 0.1")
  )

  expect_error(
    checkCriteria(criteria, project_dir = test_dir, verbose = FALSE),
    "failed validation"
  )
})

test_that("checkCriteria errors on missing column", {
  test_dir <- file.path(tempdir(), "test_check_badcol")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  supp_data <- data.table::data.table(
    gene_symbol = c("Myc"),
    value = c(0.01)
  )
  data.table::fwrite(
    supp_data,
    file.path(test_dir, ".locusPackRat", "supplementary", "test_table.csv")
  )

  criteria <- list(
    "Bad Col" = list(table = "test_table", column = "nonexistent", condition = "< 0.1")
  )

  expect_error(
    checkCriteria(criteria, project_dir = test_dir, verbose = FALSE),
    "failed validation"
  )
})

test_that("generateUpsetPlot errors on empty criteria", {
  test_dir <- file.path(tempdir(), "test_upset_empty")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  expect_error(
    generateUpsetPlot(criteria = list(), project_dir = test_dir),
    "non-empty"
  )
})

test_that("generateUpsetPlot errors on unnamed criteria", {
  test_dir <- file.path(tempdir(), "test_upset_unnamed")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  criteria <- list(
    list(table = "x", column = "y", condition = "< 0.1")
  )

  expect_error(
    generateUpsetPlot(criteria = criteria, project_dir = test_dir),
    "must be named"
  )
})

test_that("checkCriteria supports 'exists' condition", {
  test_dir <- file.path(tempdir(), "test_check_exists")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53", "Egfr"))

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  supp_data <- data.table::data.table(
    gene_symbol = c("Myc", "Tp53"),
    value = c(0.01, 0.5)
  )
  data.table::fwrite(
    supp_data,
    file.path(test_dir, ".locusPackRat", "supplementary", "test_table.csv")
  )

  criteria <- list(
    "Exists" = list(table = "test_table", column = "gene_symbol", condition = "exists")
  )

  expect_true(checkCriteria(criteria, project_dir = test_dir, verbose = FALSE))
})

test_that("checkCriteria supports 'matches' condition", {
  test_dir <- file.path(tempdir(), "test_check_matches")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53", "Egfr"))

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  supp_data <- data.table::data.table(
    gene_symbol = c("Myc", "Tp53", "Egfr"),
    category = c("cardiac", "liver", "heart disease")
  )
  data.table::fwrite(
    supp_data,
    file.path(test_dir, ".locusPackRat", "supplementary", "test_table.csv")
  )

  criteria <- list(
    "Heart Related" = list(
      table = "test_table",
      column = "category",
      condition = "matches cardiac|heart"
    )
  )

  expect_true(checkCriteria(criteria, project_dir = test_dir, verbose = FALSE))
})
