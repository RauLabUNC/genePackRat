# Tests for generateLocusZoomPlot

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

  # Test that citation message is printed
  expect_message(
    result <- generateLocusZoomPlot(
      region_id = "test_region",
      project_dir = test_dir,
      scan_table = "test_scan",
      threshold = 4,
      output_file = "minimal_zoom.pdf"
    ),
    regexp = "Plotgardener.*Bioinformatics.*2022"
  )

  expect_true(file.exists(
    file.path(test_dir, ".locusPackRat", "output", "minimal_zoom.pdf")
  ))
})
