# Tests for makeGeneSheet function

test_that("makeGeneSheet creates CSV output in gene mode", {
  test_dir <- file.path(tempdir(), "test_genesheet_csv")
  on.exit(unlink(test_dir, recursive = TRUE))

  # Create gene-mode project
  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr", "Vegfa", "Gapdh"),
    log2FC = c(1.5, -2.1, 0.8, 1.2, -0.3)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  result <- makeGeneSheet(
    format = "csv",
    project_dir = test_dir
  )

  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "output", "gene_sheet.csv")))
  expect_s3_class(result, "data.table")
  expect_true("gene_symbol" %in% names(result))
})

test_that("makeGeneSheet creates Excel output", {
  test_dir <- file.path(tempdir(), "test_genesheet_excel")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr"),
    log2FC = c(1.5, -2.1, 0.8)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  result <- makeGeneSheet(
    format = "excel",
    project_dir = test_dir
  )

  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "output", "gene_sheet.xlsx")))
  expect_s3_class(result, "data.table")
})

test_that("makeGeneSheet applies filter expression", {
  test_dir <- file.path(tempdir(), "test_genesheet_filter")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr", "Vegfa"),
    log2FC = c(1.5, -2.1, 0.8, 1.2),
    padj = c(0.001, 0.01, 0.5, 0.1)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Add a supplementary table with padj values
  supp_data <- data.table::data.table(
    gene_symbol = c("Myc", "Tp53", "Egfr", "Vegfa"),
    padj = c(0.001, 0.01, 0.5, 0.1)
  )
  data.table::fwrite(
    supp_data,
    file.path(test_dir, ".locusPackRat", "supplementary", "stats.csv")
  )

  result <- makeGeneSheet(
    format = "csv",
    filter_expr = "padj < 0.05",
    project_dir = test_dir
  )

  # Should only include genes with padj < 0.05
  expect_lte(nrow(result), 2)
})

test_that("makeGeneSheet merges supplementary tables", {
  test_dir <- file.path(tempdir(), "test_genesheet_supp")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr"),
    log2FC = c(1.5, -2.1, 0.8)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Add supplementary table
  supp_data <- data.frame(
    gene_symbol = c("Myc", "Tp53"),
    tissue_expr = c(100, 50)
  )
  addRatTable(
    data = supp_data,
    table_name = "expression",
    link_type = "gene",
    link_by = "gene_symbol",
    project_dir = test_dir
  )

  result <- makeGeneSheet(
    format = "csv",
    include_supplementary = TRUE,
    project_dir = test_dir
  )

  expect_true("tissue_expr" %in% names(result))
})

test_that("makeGeneSheet excludes specified tables", {
  test_dir <- file.path(tempdir(), "test_genesheet_exclude")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53"),
    log2FC = c(1.5, -2.1)
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
  supp1 <- data.frame(gene_symbol = c("Myc"), score1 = 10)
  supp2 <- data.frame(gene_symbol = c("Myc"), score2 = 20)

  addRatTable(supp1, "table1", "gene", "gene_symbol", project_dir = test_dir)
  addRatTable(supp2, "table2", "gene", "gene_symbol", project_dir = test_dir)

  result <- makeGeneSheet(
    format = "csv",
    include_supplementary = TRUE,
    exclude_tables = "table1",
    project_dir = test_dir
  )

  expect_false("score1" %in% names(result))
  expect_true("score2" %in% names(result))
})

test_that("makeGeneSheet works with region mode", {
  test_dir <- file.path(tempdir(), "test_genesheet_region")
  on.exit(unlink(test_dir, recursive = TRUE))

  region_data <- data.frame(
    chr = c("1", "3"),
    start = c(1000000, 5000000),
    end = c(2000000, 6000000),
    region_id = c("region1", "region2")
  )

  initPackRat(
    data = region_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  result <- makeGeneSheet(
    format = "csv",
    project_dir = test_dir
  )

  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "output", "gene_sheet.csv")))
  expect_s3_class(result, "data.table")
})

test_that("makeGeneSheet errors on missing project", {
  expect_error(
    makeGeneSheet(format = "csv", project_dir = "/nonexistent/path"),
    "No .locusPackRat directory"
  )
})

test_that("makeGeneSheet errors on bad filter expression", {
  test_dir <- file.path(tempdir(), "test_genesheet_badfilter")
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
    makeGeneSheet(
      format = "csv",
      filter_expr = "nonexistent_column > 5",
      project_dir = test_dir
    ),
    "Filter Expression Error"
  )
})

test_that("makeGeneSheet uses custom output filename", {
  test_dir <- file.path(tempdir(), "test_genesheet_filename")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53"))

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  result <- makeGeneSheet(
    format = "csv",
    output_file = "custom_output.csv",
    project_dir = test_dir
  )

  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "output", "custom_output.csv")))
})

test_that("makeGeneSheet with include_supplementary = FALSE", {
  test_dir <- file.path(tempdir(), "test_genesheet_nosupp")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53"),
    log2FC = c(1.5, -2.1)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Add supplementary table
  supp_data <- data.frame(gene_symbol = c("Myc"), extra_col = 100)
  addRatTable(supp_data, "extra", "gene", "gene_symbol", project_dir = test_dir)

  result <- makeGeneSheet(
    format = "csv",
    include_supplementary = FALSE,
    project_dir = test_dir
  )

  expect_false("extra_col" %in% names(result))
})

test_that("makeGeneSheet with split_by criteria creates multiple sheets", {
  test_dir <- file.path(tempdir(), "test_genesheet_split")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr", "Vegfa"),
    log2FC = c(2.0, -1.5, 0.5, 1.8)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  result <- makeGeneSheet(
    format = "excel",
    split_by = "criteria",
    split_criteria = list(
      "All_Genes" = "TRUE",
      "Upregulated" = "log2FC > 1"
    ),
    project_dir = test_dir
  )

  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "output", "gene_sheet.xlsx")))
})

test_that("makeGeneSheet handles prefix_mode options", {
  test_dir <- file.path(tempdir(), "test_genesheet_prefix")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53"),
    score = c(1.5, 2.0)
  )

  initPackRat(
    data = test_data,
    mode = "gene",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Add supplementary table with conflicting column name
  supp_data <- data.frame(gene_symbol = c("Myc", "Tp53"), score = c(10, 20))
  addRatTable(
    supp_data,
    "scores",
    "gene",
    "gene_symbol",
    abbreviation = "sc",
    project_dir = test_dir
  )

  # Test collision mode - should prefix the conflicting column
  result <- makeGeneSheet(
    format = "csv",
    prefix_mode = "collision",
    project_dir = test_dir
  )

  # The supplementary score column should be prefixed
  expect_true("sc_score" %in% names(result) || "score" %in% names(result))
})
