# Extended tests for addRatTable function - point and region link types

test_that("addRatTable handles point data (link_type = 'point')", {
  test_dir <- file.path(tempdir(), "test_add_point")
  on.exit(unlink(test_dir, recursive = TRUE))

  # Create region-mode project
  region_data <- data.frame(
    chr = c("1", "1"),
    start = c(1000000, 2000000),
    end = c(1500000, 2500000),
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

  # Add SNP-style point data
  snp_data <- data.frame(
    chr = c("1", "1", "1"),
    pos = c(1200000, 1400000, 2200000),
    snp_id = c("rs123", "rs456", "rs789"),
    pvalue = c(0.001, 0.05, 0.01)
  )

  result <- addRatTable(
    data = snp_data,
    table_name = "snps",
    link_type = "point",
    project_dir = test_dir
  )

  expect_true(result)
  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "supplementary", "snps.csv")))

  # Verify the saved data has start/end columns
  saved_data <- data.table::fread(
    file.path(test_dir, ".locusPackRat", "supplementary", "snps.csv")
  )
  expect_true("start" %in% names(saved_data))
  expect_true("end" %in% names(saved_data))
})

test_that("addRatTable point data normalizes position columns", {
  test_dir <- file.path(tempdir(), "test_add_point_normalize")
  on.exit(unlink(test_dir, recursive = TRUE))

  region_data <- data.frame(
    chr = "1",
    start = 1000000,
    end = 2000000
  )

  initPackRat(
    data = region_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Use 'bp' instead of 'pos'
  snp_data <- data.frame(
    chrom = "1",
    bp = 1500000,
    snp_id = "rs999"
  )

  result <- addRatTable(
    data = snp_data,
    table_name = "snps_bp",
    link_type = "point",
    project_dir = test_dir
  )

  expect_true(result)

  saved_data <- data.table::fread(
    file.path(test_dir, ".locusPackRat", "supplementary", "snps_bp.csv")
  )
  # Should have chr (normalized from chrom) and pos (normalized from bp)
  expect_true("chr" %in% names(saved_data))
  expect_true("pos" %in% names(saved_data))
})

test_that("addRatTable handles region data (link_type = 'region')", {
  test_dir <- file.path(tempdir(), "test_add_region")
  on.exit(unlink(test_dir, recursive = TRUE))

  # Create region-mode project
  region_data <- data.frame(
    chr = c("1", "3"),
    start = c(1000000, 5000000),
    end = c(2000000, 6000000),
    region_id = c("peak1", "peak2")
  )

  initPackRat(
    data = region_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Add region-linked supplementary data
  supp_regions <- data.frame(
    chr = c("1", "3"),
    start = c(1000000, 5000000),
    end = c(2000000, 6000000),
    score = c(10, 20)
  )

  result <- addRatTable(
    data = supp_regions,
    table_name = "scores",
    link_type = "region",
    link_by = "chr,start,end",
    project_dir = test_dir
  )

  expect_true(result)
  expect_true(file.exists(file.path(test_dir, ".locusPackRat", "supplementary", "scores.csv")))
})

test_that("addRatTable with region mode in gene project extracts genes", {
  test_dir <- file.path(tempdir(), "test_add_region_gene")
  on.exit(unlink(test_dir, recursive = TRUE))

  # Create gene-mode project
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

  # Add gene-linked data
  gene_data <- data.frame(
    gene_symbol = c("Myc", "Tp53"),
    tissue = c("liver", "brain")
  )

  result <- addRatTable(
    data = gene_data,
    table_name = "tissues",
    link_type = "gene",
    link_by = "gene_symbol",
    project_dir = test_dir
  )

  expect_true(result)

  # Check config was updated
  config <- jsonlite::read_json(file.path(test_dir, ".locusPackRat", "config.json"))
  expect_equal(config$supplementary_tables$tissues$link_type, "gene")
})

test_that("addRatTable point data errors without position column", {
  test_dir <- file.path(tempdir(), "test_add_point_err")
  on.exit(unlink(test_dir, recursive = TRUE))

  region_data <- data.frame(
    chr = "1",
    start = 1000000,
    end = 2000000
  )

  initPackRat(
    data = region_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Missing position column
  bad_data <- data.frame(
    chr = "1",
    snp_id = "rs999"
  )

  expect_error(
    addRatTable(
      data = bad_data,
      table_name = "bad_snps",
      link_type = "point",
      project_dir = test_dir
    ),
    "position column"
  )
})

test_that("addRatTable point data errors without chromosome column", {
  test_dir <- file.path(tempdir(), "test_add_point_nochr")
  on.exit(unlink(test_dir, recursive = TRUE))

  region_data <- data.frame(
    chr = "1",
    start = 1000000,
    end = 2000000
  )

  initPackRat(
    data = region_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Missing chromosome column
  bad_data <- data.frame(
    pos = 1500000,
    snp_id = "rs999"
  )

  expect_error(
    addRatTable(
      data = bad_data,
      table_name = "bad_snps",
      link_type = "point",
      project_dir = test_dir
    ),
    "chromosome column"
  )
})

test_that("addRatTable region mode links to existing regions", {
  test_dir <- file.path(tempdir(), "test_add_region_link")
  on.exit(unlink(test_dir, recursive = TRUE))

  # Create region-mode project
  region_data <- data.frame(
    chr = c("1"),
    start = c(1000000),
    end = c(2000000),
    region_id = c("test_region")
  )

  initPackRat(
    data = region_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  # Add overlapping region data
  overlap_data <- data.frame(
    chr = c("1"),
    start = c(1500000),
    end = c(1800000),
    feature = c("enhancer")
  )

  result <- addRatTable(
    data = overlap_data,
    table_name = "features",
    link_type = "region",
    link_by = "chr,start,end",
    project_dir = test_dir
  )

  expect_true(result)

  # Check metadata
  config <- jsonlite::read_json(file.path(test_dir, ".locusPackRat", "config.json"))
  expect_equal(config$supplementary_tables$features$link_type, "region")
  expect_equal(config$supplementary_tables$features$link_by, "chr,start,end")
})

test_that("addRatTable stores original link_type for point data", {
  test_dir <- file.path(tempdir(), "test_add_point_config")
  on.exit(unlink(test_dir, recursive = TRUE))

  region_data <- data.frame(
    chr = "1",
    start = 1000000,
    end = 2000000
  )

  initPackRat(
    data = region_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    force = TRUE
  )

  snp_data <- data.frame(
    chr = "1",
    pos = 1500000,
    snp_id = "rs123"
  )

  addRatTable(
    data = snp_data,
    table_name = "test_snps",
    link_type = "point",
    project_dir = test_dir
  )

  # Config should store "point" as the original link type
  config <- jsonlite::read_json(file.path(test_dir, ".locusPackRat", "config.json"))
  expect_equal(config$supplementary_tables$test_snps$link_type, "point")
})
