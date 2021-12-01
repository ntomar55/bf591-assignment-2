#!/usr/bin/Rscript
# if you chnge the name of your script, this line must be changed as well
source("main.R")
library(testthat)

test_that("loading csv correctly", {
  result_tib <- load_expression("/project/bf528/project_1/data/example_intensity_data.csv")
  expect_equal(dim(result_tib), c(54675, 36))
  expect_true(is_tibble(result_tib))
})

test_that("plot_ggplot() correctly creates a boxplot from sample data", {
  plot_tib <- tibble(probeids = c("202274_at", "202541_at", "202542_s_at", "203919_at"),
                     hgnc = c("ACTG2", "AIMP1", "AIMP1", "TCEA2"),
                     gene_set = rep("good", 4),
                     GSM1 = c(8.05, 8.40, 9.55, 4.44),
                     GSM2 = c(7.74, 7.11, 8.48, 5.39))
  p <- plot_ggplot(plot_tib)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomBoxplot")
})