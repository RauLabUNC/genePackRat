# Tests for listPackRatTables and removeRatTable functions

test_that("listPackRatTables returns empty for new project", {
  test_dir <- file.path(tempdir(), "test_list_empty")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  result <- listPackRatTables(project_dir = test_dir)

  expect_equal(nrow(result), 0)
})

test_that("listPackRatTables lists added tables", {
  test_dir <- file.path(tempdir(), "test_list_tables")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  supp_data <- data.frame(gene_symbol = c("Myc"), score = 10)
  addRatTable(supp_data, "scores", "gene", "gene_symbol", project_dir = test_dir)

  result <- listPackRatTables(project_dir = test_dir)

  expect_equal(nrow(result), 1)
  expect_equal(result$table_name[1], "scores")
})

test_that("listPackRatTables requires existing project", {
  expect_error(
    listPackRatTables(project_dir = "/nonexistent/path"),
    "No .locusPackRat directory"
  )
})

test_that("removeRatTable removes table and config entry", {
  test_dir <- file.path(tempdir(), "test_remove")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  supp_data <- data.frame(gene_symbol = c("Myc"), score = 10)
  addRatTable(supp_data, "to_remove", "gene", "gene_symbol", project_dir = test_dir)

  # Verify table exists
  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "supplementary", "to_remove.csv")))

  # Remove it
  removeRatTable("to_remove", project_dir = test_dir)

  # Verify table is gone
  expect_false(file.exists(file.path(test_dir, ".locusPackRat", "supplementary", "to_remove.csv")))

  config <- jsonlite::read_json(file.path(test_dir, ".locusPackRat", "config.json"))
  expect_null(config$supplementary_tables$to_remove)
})

test_that("removeRatTable errors for nonexistent table", {
  test_dir <- file.path(tempdir(), "test_remove_err")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  expect_error(
    removeRatTable("nonexistent", project_dir = test_dir),
    "not found"
  )
})

test_that("listPackRatTables with full_info prints column names", {
  test_dir <- file.path(tempdir(), "test_list_fullinfo")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  supp_data <- data.frame(
    gene_symbol = c("Myc", "Tp53"),
    score = c(10, 20),
    category = c("A", "B")
  )
  addRatTable(supp_data, "test_table", "gene", "gene_symbol", project_dir = test_dir)

  # Should run without error and return table info
  expect_message(
    result <- listPackRatTables(project_dir = test_dir, full_info = TRUE),
    "Printing full column names"
  )

  expect_equal(nrow(result), 1)
})

test_that("listPackRatTables with criteria_info prints column summaries", {
  test_dir <- file.path(tempdir(), "test_list_criteria")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53", "Egfr"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  supp_data <- data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr"),
    pvalue = c(0.001, 0.05, 0.5),
    category = c("sig", "sig", "ns"),
    is_validated = c(TRUE, FALSE, TRUE)
  )
  addRatTable(supp_data, "stats", "gene", "gene_symbol", project_dir = test_dir)

  # Should print detailed column info
  expect_message(
    result <- listPackRatTables(project_dir = test_dir, criteria_info = TRUE),
    "Column Details for Criteria Building"
  )

  expect_equal(nrow(result), 1)
})

test_that("listPackRatTables criteria_info shows numeric ranges", {
  test_dir <- file.path(tempdir(), "test_list_criteria_num")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  supp_data <- data.frame(
    gene_symbol = c("Myc", "Tp53"),
    score = c(10.5, 20.3)
  )
  addRatTable(supp_data, "scores", "gene", "gene_symbol", project_dir = test_dir)

  # Should show min, max, median for numeric columns
  expect_message(
    listPackRatTables(project_dir = test_dir, criteria_info = TRUE),
    "numeric"
  )
})

test_that("listPackRatTables criteria_info shows character unique values", {
  test_dir <- file.path(tempdir(), "test_list_criteria_char")
  on.exit(unlink(test_dir, recursive = TRUE))

  test_data <- data.frame(gene_symbol = c("Myc", "Tp53", "Egfr"))
  initPackRat(test_data, "gene", "mouse", "mm39", test_dir, force = TRUE)

  supp_data <- data.frame(
    gene_symbol = c("Myc", "Tp53", "Egfr"),
    tissue = c("liver", "brain", "liver")
  )
  addRatTable(supp_data, "tissues", "gene", "gene_symbol", project_dir = test_dir)

  # Should show unique values for character columns
  expect_message(
    listPackRatTables(project_dir = test_dir, criteria_info = TRUE),
    "character"
  )
})
