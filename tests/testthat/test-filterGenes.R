# Tests for filterGenes and makeFilter functions

test_that("makeFilter creates valid filter specification", {
  f <- makeFilter("score", ">", 50)

  expect_equal(f$column, "score")
  expect_equal(f$condition, ">")
  expect_equal(f$value, 50)
})

test_that("makeFilter validates condition", {
  expect_error(makeFilter("col", "invalid", 1), "Invalid condition")
})

test_that("filterGenes with semi join keeps matching rows", {
  main <- data.frame(
    gene = c("A", "B", "C", "D"),
    value = c(1, 2, 3, 4)
  )

  filter_df <- data.frame(gene = c("A", "C"))

  result <- filterGenes(main, filter_df, joinType = "semi", by = "gene")

  expect_equal(nrow(result), 2)
  expect_true(all(result$gene %in% c("A", "C")))
})

test_that("filterGenes with anti join removes matching rows", {
  main <- data.frame(
    gene = c("A", "B", "C", "D"),
    value = c(1, 2, 3, 4)
  )

  filter_df <- data.frame(gene = c("A", "C"))

  result <- filterGenes(main, filter_df, joinType = "anti", by = "gene")

  expect_equal(nrow(result), 2)
  expect_true(all(result$gene %in% c("B", "D")))
})

test_that("filterGenes applies filter conditions", {
  main <- data.frame(
    gene = c("A", "B", "C", "D"),
    score = c(10, 50, 30, 80)
  )

  filter_df <- data.frame(gene = c("A", "B", "C", "D"))
  filters <- list(makeFilter("score", ">", 25))

  result <- filterGenes(main, filter_df, joinType = "inner", by = "gene", filters = filters)

  expect_true(all(result$score > 25))
})

test_that("filterGenes handles empty result", {
  main <- data.frame(
    gene = c("A", "B"),
    value = c(1, 2)
  )

  filter_df <- data.frame(gene = c("X", "Y"))

  result <- filterGenes(main, filter_df, joinType = "semi", by = "gene")

  expect_equal(nrow(result), 0)
})

test_that("filterGenes with %in% condition", {
  main <- data.frame(
    gene = c("A", "B", "C", "D"),
    category = c("good", "bad", "good", "ugly")
  )

  filter_df <- data.frame(gene = c("A", "B", "C", "D"))
  filters <- list(makeFilter("category", "%in%", c("good", "ugly")))

  result <- filterGenes(main, filter_df, joinType = "inner", by = "gene", filters = filters)

  expect_equal(nrow(result), 3)
  expect_true(all(result$category %in% c("good", "ugly")))
})
