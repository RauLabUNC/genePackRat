# Testing Plan for Untested Functions

## Overview
Two main functions remain untested:
1. `generateGeneInfoExcel()` - Creates 3-sheet Excel workbook (629-1018 lines)
2. `generateLocusPacket()` - Orchestrates full packet generation (1045-1300+ lines)

## Function Analysis

### generateGeneInfoExcel()
**Purpose**: Creates Excel workbook with gene annotations across 3 sheets

**Inputs Required**:
- genes_in_locus (data.table) - Genes in region
- loci_info (data.frame) - QTL locus information
- merged_gene_info (data.table) - Pre-merged gene summaries
- genes_mouse (data.table) - Gene coordinates
- ortho_mouse2h (data.table) - Mouse-human orthologs
- mouse_pheno (data.table) - MGI phenotypes
- associations (data.table) - Human disease associations
- rna_info (data.table) - Expression data
- founder_mutations (data.table, optional) - CC variants
- output_file (string) - Path for Excel file

**Operations**:
1. Filters genes to locus
2. Merges multiple data sources
3. Calculates mean expression by condition
4. Creates T/F columns for each locus
5. Generates 3 Excel sheets with formatting

### generateLocusPacket()
**Purpose**: Complete packet generation orchestrator

**Inputs Required**:
- locus_cluster (data.frame) - Locus definitions
- scan_data, threshold_data - QTL results
- 9+ file paths for data sources
- output_path - Directory for results

**Operations**:
1. Loads 7-9 CSV files
2. Creates output directory structure
3. Calls generateGeneInfoExcel()
4. Calls generateLocusZoomPlot() for each QTL
5. Creates README summary
6. Optionally creates SNP table

## Testing Strategy

### Phase 1: Minimal Mock Data Creation

#### Create 10-Gene Test Dataset
```r
create_test_data_for_excel <- function() {
  # 10 genes on chr 1 from 100-110 Mb
  genes_mouse <- data.table(
    mouse_ensembl_id = paste0("ENSMUSG0000000", 1:10),
    mouse_gene_symbol = paste0("TestGene", 1:10),
    chr = 1,
    start_bp = seq(100e6, 109e6, 1e6),
    end_bp = seq(100.5e6, 109.5e6, 1e6),
    gene_biotype = "protein_coding"
  )

  # Merged gene info with required columns
  merged_gene_info <- data.table(
    mouse_ensembl_id = genes_mouse$mouse_ensembl_id,
    mouse_gene_symbol = genes_mouse$mouse_gene_symbol,
    n_trait_drug = sample(0:5, 10, replace=TRUE),
    trait_drug = rep("HR_Ctrl,BW_Iso", 10),
    avgenes_cpm = runif(10, 0, 100),
    disease = rep("cardiovascular disease", 10),
    ontology = rep("abnormal heart morphology", 10),
    pubmedID = rep("12345678", 10)
  )

  # Expression data (samples x genes)
  rna_info <- data.table(
    Sample = paste0("S", 1:8),
    Drug = rep(c("Ctrl", "Iso"), 4),
    Sex = rep(c("M", "F"), each=4),
    TestGene1 = rnorm(8, 50, 10),
    TestGene2 = rnorm(8, 30, 5),
    # ... more genes
  )

  # Mouse phenotypes
  mouse_pheno <- data.table(
    mouse_gene_symbol = rep(genes_mouse$mouse_gene_symbol[1:5], each=2),
    OntologyAnnotation.ontologyTerm.identifier = paste0("MP:000", 1:10),
    OntologyAnnotation.ontologyTerm.name = rep(c("heart", "vessel"), 5),
    OntologyAnnotation.evidence.publications.pubMedId = "12345678",
    OntologyAnnotation.evidence.comments.description = "test phenotype"
  )

  # Human orthologs
  ortho_mouse2h <- data.table(
    mouse_ensembl_id = genes_mouse$mouse_ensembl_id,
    human_ensembl_id = paste0("ENSG0000000", 1:10),
    human_gene_symbol = paste0("HUMAN", 1:10)
  )

  # Disease associations
  associations <- data.table(
    human_ensembl_id = ortho_mouse2h$human_ensembl_id,
    symbol = ortho_mouse2h$human_gene_symbol,
    disease_id = "MONDO_0004995",
    disease_name = "cardiovascular disease",
    association_score = runif(10, 0, 1)
  )

  list(
    genes_mouse = genes_mouse,
    merged_gene_info = merged_gene_info,
    rna_info = rna_info,
    mouse_pheno = mouse_pheno,
    ortho_mouse2h = ortho_mouse2h,
    associations = associations
  )
}
```

### Phase 2: Unit Tests for generateGeneInfoExcel()

#### Test 1: Basic Excel Generation
```r
test_that("generateGeneInfoExcel creates Excel with correct structure", {
  # Create test data
  test_data <- create_test_data_for_excel()

  # Define test locus
  loci_info <- data.frame(
    chr = 1,
    upper_pos_lod_drop = 100,
    lower_pos_lod_drop = 110,
    peak_pos = 105,
    trait = "HR",
    drug = "Ctrl"
  )

  # Genes in locus (all 10 test genes)
  genes_in_locus <- test_data$genes_mouse

  # Create temp file
  temp_excel <- tempfile(fileext = ".xlsx")

  # Generate Excel
  result <- generateGeneInfoExcel(
    genes_in_locus = genes_in_locus,
    loci_info = loci_info,
    merged_gene_info = test_data$merged_gene_info,
    genes_mouse = test_data$genes_mouse,
    ortho_mouse2h = test_data$ortho_mouse2h,
    mouse_pheno = test_data$mouse_pheno,
    associations = test_data$associations,
    rna_info = test_data$rna_info,
    founder_mutations = NULL,
    output_file = temp_excel
  )

  # Verify file exists
  expect_true(file.exists(temp_excel))

  # Read back and verify structure
  wb <- openxlsx::loadWorkbook(temp_excel)
  sheets <- openxlsx::sheets(wb)

  expect_equal(length(sheets), 3)
  expect_true("AllGenesInCluster" %in% sheets)
  expect_true("AllMousePhenotypes" %in% sheets)
  expect_true("AllHumanDiseases" %in% sheets)

  # Check main sheet content
  main_sheet <- openxlsx::readWorkbook(wb, sheet = 1)
  expect_equal(nrow(main_sheet), 10)  # 10 test genes
  expect_true("Mouse Gene Symbol" %in% names(main_sheet))
  expect_true("chr" %in% names(main_sheet))

  # Cleanup
  unlink(temp_excel)
})
```

#### Test 2: Empty Locus Handling
```r
test_that("generateGeneInfoExcel handles empty locus gracefully", {
  test_data <- create_test_data_for_excel()

  # Empty genes_in_locus
  genes_in_locus <- test_data$genes_mouse[0, ]

  result <- generateGeneInfoExcel(
    genes_in_locus = genes_in_locus,
    # ... other params
    output_file = tempfile(fileext = ".xlsx")
  )

  expect_null(result)
})
```

#### Test 3: Expression Calculations
```r
test_that("generateGeneInfoExcel calculates expression means correctly", {
  # Test that Drug_Sex combinations are calculated properly
  # Verify Ave_Exp_Ctrl_F, Ave_Exp_Ctrl_M, etc. columns
})
```

### Phase 3: Integration Tests for generateLocusPacket()

#### Test with Mocked File I/O
```r
test_that("generateLocusPacket orchestrates correctly with mock data", {
  # Use mockery or testthat mocking
  with_mocked_bindings(
    fread = function(file) {
      if (grepl("genes_mouse", file)) return(test_genes_mouse)
      if (grepl("orthology", file)) return(test_ortho)
      # ... etc
    },
    code = {
      locus_cluster <- data.frame(
        chr = 1,
        upper_pos_lod_drop = 100,
        lower_pos_lod_drop = 110,
        peak_pos = 105,
        trait = "HR",
        drug = "Ctrl"
      )

      # Mock scan_data and threshold_data
      scan_data <- list(HR_Ctrl = test_scan_data)
      threshold_data <- list(HR_Ctrl_threshold = 5.0)

      output_dir <- tempdir()

      result <- generateLocusPacket(
        locus_cluster = locus_cluster,
        input_path = "fake/path",
        output_path = output_dir,
        scan_data = scan_data,
        threshold_data = threshold_data
      )

      # Verify directory structure created
      expect_true(dir.exists(file.path(output_dir, "locus_chr1_100-110Mb")))
    }
  )
})
```

### Phase 4: Refactoring Recommendations

#### Extract Business Logic from I/O
```r
# Separate calculation from Excel writing
calculate_gene_summary <- function(genes_in_locus, merged_gene_info, ...) {
  # Pure data transformation, no I/O
  # Returns list of data.tables for sheets
}

write_excel_sheets <- function(sheet_data, output_file) {
  # Pure I/O operation
  # Takes calculated data and writes Excel
}
```

#### Make Directory Operations Testable
```r
# Extract directory creation
create_packet_structure <- function(output_dir, locus_name) {
  dirs <- list(
    main = file.path(output_dir, locus_name),
    plots = file.path(output_dir, locus_name, "zoomPlots")
  )
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  return(dirs)
}
```

## Implementation Priority

### Week 1: Core Testing
1. ✅ Create minimal 10-gene test dataset
2. ✅ Write 3 basic tests for generateGeneInfoExcel
3. ✅ Verify Excel structure tests pass

### Week 2: Advanced Testing
4. Add expression calculation tests
5. Add locus presence (T/F) column tests
6. Test error handling paths

### Week 3: Integration
7. Mock file I/O for generateLocusPacket
8. Test directory creation
9. Test function orchestration

## Testing Challenges & Solutions

### Challenge 1: Excel Binary Format
**Solution**: Test structure, not bytes
- Verify sheet names exist
- Check row/column counts
- Sample specific cells

### Challenge 2: Large Dependencies
**Solution**: Create minimal valid dataset
- 10 genes instead of 2000+
- Maintain all relationships
- ~1MB instead of 65MB

### Challenge 3: File I/O Coupling
**Solution**: Mock at boundaries
- Mock fread() for input
- Use tempfile() for output
- Verify file existence, not content

### Challenge 4: Complex Data Relationships
**Solution**: Test fixtures with known results
- Pre-calculate expected means
- Use simple test data (1,2,3...)
- Document expected transformations

## Success Criteria

✅ Tests run in < 5 seconds
✅ No real file dependencies
✅ 80% code coverage for logic
✅ All error paths tested
✅ Works on CI/CD systems

## Next Steps

1. Create `tests/testthat/fixtures/excel_test_data.R`
2. Implement first 3 tests
3. Run and iterate
4. Add to CI pipeline
5. Document test data schema