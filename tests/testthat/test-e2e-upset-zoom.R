# End-to-end test for generateUpsetPlot and generateLocusZoomPlot
# Based on scripts/demo_upset_plot.R but with mock data (no API calls)

test_that("End-to-end upset and locus zoom workflow completes successfully", {
  skip_on_cran()

  test_dir <- file.path(tempdir(), "test_e2e_upset_zoom")
  on.exit(unlink(test_dir, recursive = TRUE))

  # Initialize project with chr3 region
  region_data <- data.table::data.table(
    chr = "3",
    start = 131000000,
    end = 150000000,
    region_id = "chr3_locus"
  )

  initPackRat(
    data = region_data,
    mode = "region",
    species = "mouse",
    genome = "mm39",
    project_dir = test_dir,
    keep_pseudo = FALSE,
    overlap_mode = "any",
    force = TRUE
  )

  expect_true(dir.exists(file.path(test_dir, ".locusPackRat")))

  # Get detected genes from the initialized project
  regions_dt <- data.table::fread(
    file.path(test_dir, ".locusPackRat", "input", "regions.csv")
  )
  detected_genes <- unlist(strsplit(regions_dt$genes, ",\\s*"))
  n_genes <- length(detected_genes)

  expect_gt(n_genes, 0)

  # --- Add simulated differential expression results ---
  set.seed(123)
  deseq_results <- data.table::data.table(
    gene_symbol = detected_genes,
    log2FoldChange = rnorm(n_genes, mean = 0, sd = 1.5),
    pvalue = runif(n_genes, 0, 1),
    padj = NA_real_,
    baseMean = abs(rnorm(n_genes, 500, 300))
  )
  # Make ~30% of genes significantly DE
  sig_idx <- sample(seq_len(n_genes), round(n_genes * 0.3))
  deseq_results[sig_idx, pvalue := runif(.N, 0, 0.01)]
  deseq_results[, padj := p.adjust(pvalue, method = "BH")]

  addRatTable(
    data = deseq_results,
    table_name = "deseq_results",
    abbreviation = "de",
    link_type = "gene",
    project_dir = test_dir
  )

  # --- Simulated eQTL data (~40% of genes) ---
  eqtl_genes <- sample(detected_genes, round(n_genes * 0.4))
  eqtl_data <- data.table::data.table(
    gene_symbol = eqtl_genes,
    snp_id = paste0("rs", sample(100000:999999, length(eqtl_genes))),
    pvalue = runif(length(eqtl_genes), 0, 0.05),
    effect_size = rnorm(length(eqtl_genes), 0, 0.5)
  )

  addRatTable(
    data = eqtl_data,
    table_name = "eqtl_data",
    abbreviation = "eq",
    link_type = "gene",
    project_dir = test_dir
  )

  # --- Simulated CC variant data (~25% of genes) ---
  cc_genes <- sample(detected_genes, round(n_genes * 0.25))
  cc_variants <- data.table::data.table(
    gene_symbol = cc_genes,
    variant_type = sample(c("missense", "nonsense", "frameshift"),
                          length(cc_genes), replace = TRUE),
    is_deleterious = sample(c(TRUE, FALSE), length(cc_genes),
                            replace = TRUE, prob = c(0.7, 0.3))
  )

  addRatTable(
    data = cc_variants,
    table_name = "cc_variants",
    abbreviation = "cc",
    link_type = "gene",
    project_dir = test_dir
  )

  # --- Mock mouse_phenotypes.csv (instead of queryMouseMine) ---
  phenotype_genes <- sample(detected_genes, round(n_genes * 0.4))
  mock_phenotypes <- data.table::data.table(
    gene_symbol = phenotype_genes,
    phenotype = sample(
      c("cardiac hypertrophy", "heart rate variability", "ventricular defect",
        "liver abnormality", "kidney dysfunction", "neurological phenotype"),
      length(phenotype_genes),
      replace = TRUE
    ),
    mp_id = paste0("MP:", sprintf("%07d", sample(1:99999, length(phenotype_genes))))
  )
  data.table::fwrite(
    mock_phenotypes,
    file.path(test_dir, ".locusPackRat", "supplementary", "mouse_phenotypes.csv")
  )

  # --- Add simulated scan data for locus zoom ---
  sim_positions <- seq(131000000, 150000000, by = 50000)
  scan_data <- data.table::data.table(
    chr = 3,
    pos = sim_positions,
    start = sim_positions,
    end = sim_positions,
    lod = abs(rnorm(length(sim_positions), 2, 1.5)),
    p = 10^(-abs(rnorm(length(sim_positions), 2, 1.5)))
  )
  # Add a peak
  peak_idx <- which.min(abs(sim_positions - 140000000))
  scan_data[peak_idx, lod := 6.5]

  addRatTable(
    data = scan_data,
    table_name = "scan_data",
    abbreviation = "sc",
    link_type = "region",
    link_by = "chr,start,end",
    project_dir = test_dir
  )

  # --- Test listPackRatTables ---
  tables <- listPackRatTables(test_dir, criteria_info = TRUE)
  expect_true("deseq_results" %in% tables$table_name)
  expect_true("eqtl_data" %in% tables$table_name)
  expect_true("mouse_phenotypes" %in% tables$table_name)

  # --- Define criteria for upset plot ---
  criteria <- list(
    "DE (FDR < 0.1)" = list(
      table = "deseq_results",
      column = "padj",
      condition = "< 0.1"
    ),
    "Has cis-eQTL" = list(
      table = "eqtl_data",
      column = "gene_symbol",
      condition = "exists"
    ),
    "CC Variant" = list(
      table = "cc_variants",
      column = "gene_symbol",
      condition = "exists"
    ),
    "Cardiac Phenotype (Mouse)" = list(
      table = "mouse_phenotypes",
      column = "phenotype",
      condition = "matches cardiac|heart|ventric"
    )
  )

  # --- Test checkCriteria ---
  expect_true(checkCriteria(criteria, project_dir = test_dir, verbose = FALSE))

  # --- Find genes for highlighting ---
  de_genes <- deseq_results[padj < 0.1, gene_symbol]
  cardiac_genes <- mock_phenotypes[
    grepl("cardiac|heart|ventric", phenotype, ignore.case = TRUE),
    unique(gene_symbol)
  ]
  highlight_candidates <- intersect(de_genes, cardiac_genes)

  # --- Test generateUpsetPlot ---
  upset_result <- generateUpsetPlot(
    criteria = criteria,
    project_dir = test_dir,
    highlight_genes = if (length(highlight_candidates) > 0) highlight_candidates else NULL,
    highlight_intersections = list(
      c("DE (FDR < 0.1)", "Cardiac Phenotype (Mouse)")
    ),
    highlight_color = "#a03b60",
    intersection_color = "#53add2",
    width = 9,
    height = 5.5,
    output_file = "test_upset.pdf"
  )

  expect_true(file.exists(
    file.path(test_dir, ".locusPackRat", "output", "test_upset.pdf")
  ))

  # --- Test generateLocusZoomPlot ---
  zoom_result <- generateLocusZoomPlot(
    region_id = "chr3_locus",
    project_dir = test_dir,
    scan_table = "scan_data",
    width = 10,
    height = 6,
    threshold = 4,
    highlight_genes = if (length(highlight_candidates) > 0) highlight_candidates else NULL,
    output_file = "test_locus_zoom.pdf"
  )

  expect_true(file.exists(
    file.path(test_dir, ".locusPackRat", "output", "test_locus_zoom.pdf")
  ))
})
