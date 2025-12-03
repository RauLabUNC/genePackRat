# Tests for addRatTable function

test_that("addRatTable adds gene-linked table", {
  test_dir <- file.path(tempdir(), "test_add_gene")
  on.exit(unlink(test_dir, recursive = TRUE))

  # Create project
  test_data <- data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr"),
    log2FC = c(1.5, -2.1, 0.8)
  )
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  # Add supplementary table
  supp_data <- data.frame(
    gene_symbol = c("Myc", "Tp53"),
    tissue_expr = c(100, 50)
  )

  result <- addRatTable(
    data = supp_data,
    table_name = "expression",
    link_type = "gene",
    link_by = "gene_symbol",
    project_dir = test_dir
  )

  expect_true(result)
  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "supplementary", "expression.csv")))
})

test_that("addRatTable validates table_name characters", {
  test_dir <- file.path(tempdir(), "test_add_validate")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  supp_data <- data.frame(gene_symbol = c("Myc"), value = 1)

  expect_error(
    addRatTable(supp_data, "bad@name!", "gene", "gene_symbol", project_dir = test_dir),
    "must contain only letters"
  )
})

test_that("addRatTable prevents duplicate abbreviations", {
  test_dir <- file.path(tempdir(), "test_add_abbrev")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  supp1 <- data.frame(gene_symbol = c("Myc"), value = 1)
  supp2 <- data.frame(gene_symbol = c("Myc"), score = 2)

  addRatTable(supp1, "table1", "gene", "gene_symbol", abbreviation = "ab", project_dir = test_dir)

  expect_error(
    addRatTable(supp2, "table2", "gene", "gene_symbol", abbreviation = "ab", project_dir = test_dir),
    "already used"
  )
})

test_that("addRatTable requires existing project", {
  test_dir <- file.path(tempdir(), "nonexistent_project")
  supp_data <- data.frame(gene_symbol = c("Myc"), value = 1)

  expect_error(
    addRatTable(supp_data, "test", "gene", "gene_symbol", project_dir = test_dir),
    "No .locusPackRat directory"
  )
})

test_that("addRatTable updates config with table metadata", {
  test_dir <- file.path(tempdir(), "test_add_config")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  supp_data <- data.frame(
    gene_symbol = c("Myc", "Tp53"),
    score = c(10, 20)
  )
  addRatTable(supp_data, "scores", "gene", "gene_symbol",
              abbreviation = "sc", project_dir = test_dir)

  config <- jsonlite::read_json(file.path(test_dir, ".locusPackRat", "config.json"))

  expect_true("scores" %in% names(config$supplementary_tables))
  expect_equal(config$supplementary_tables$scores$abbreviation, "sc")
  expect_equal(config$supplementary_tables$scores$link_type, "gene")
})
