test_that("plotDysbiosis works", {
  library(dysbiosisR)
  data("WirbelJ_2018")
  dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
  # get reference samples
  ref.samples <- sample_names(subset_samples(WirbelJ_2018,
                                             disease == "healthy"))

  dysbiosis_1 <- dysbiosisMedianCLV(WirbelJ_2018,
                                    dist_mat = dist.mat,
                                    reference_samples = ref.samples)
  expect_error(plotDysbiosis(df=NULL,
                             xvar="disease",
                             yvar="score",
                             colors=c(CRC="brown3", healthy="steelblue")))

  expect_error(plotDysbiosis(df=dysbiosis_1,
                             xvar="diseaseXX",
                             yvar="score",
                             colors=c(CRC="brown3", healthy="steelblue")))

  expect_error(plotDysbiosis(df=dysbiosis_1,
                             xvar="disease",
                             yvar="score.obb",
                             colors=c(CRC="brown3", healthy="steelblue")))
})
