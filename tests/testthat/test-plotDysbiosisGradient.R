test_that("plotDysbiosisGradient works", {

  library(dysbiosisR)
  volcano <- c("#003f5c", "#58508d","#bc5090","#ff6361", "#ffa600")
  dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
  # get reference samples
  ref.samples <- sample_names(subset_samples(WirbelJ_2018,
                                             disease == "healthy"))

  dysbiosis_1 <- dysbiosisMedianCLV(WirbelJ_2018,
                                    dist_mat = dist.mat,
                                    reference_samples = ref.samples)

  expect_error(plotDysbiosisGradient(df=dysbiosis_1,
                                     score="oob.score",
                                     high_line = 0.5,
                                     # low_line = normobiosis_thres,
                                     group_var = "disease",
                                     group_colors=c("healthy" = "steelblue",
                                                    "CRC"= "brown3"),
                                     point_size = 2,
                                     bg_colors = rev(volcano),
                                     jitter_width = 0.1))

  expect_error(plotDysbiosisGradient(df=NULL,
                                     score=NULL,
                                     high_line = 0.5,
                                     # low_line = normobiosis_thres,
                                     group_var = "disease",
                                     group_colors=NULL,
                                     point_size = 2,
                                     bg_colors = rev(volcano),
                                     jitter_width = 0.1))

  expect_error(plotDysbiosisGradient(df=dysbiosis_1,
                                     score="score",
                                     high_line = 0.5,
                                     # low_line = normobiosis_thres,
                                     group_var = "notfound",
                                     group_colors=c("healthy" = "steelblue",
                                                    "CRC"= "brown3"),
                                     point_size = 2,
                                     bg_colors = rev(volcano),
                                     jitter_width = 0.1))

  expect_error(plotDysbiosisGradient(df=dysbiosis_1,
                                     score="score",
                                     high_line = 0.5,
                                     # low_line = normobiosis_thres,
                                     group_var = "notfound",
                                     group_colors=NULL,
                                     point_size = 2,
                                     bg_colors = rev(volcano),
                                     jitter_width = 0.1))

  expect_error(plotDysbiosisGradient(df=NULL,
                                     score="score",
                                     high_line = 0.5,
                                     # low_line = normobiosis_thres,
                                     group_var = "disease",
                                     group_colors=NULL,
                                     point_size = 2,
                                     bg_colors = rev(volcano),
                                     jitter_width = 0.1))

  p <- plotDysbiosisGradient(df=dysbiosis_1,
                             score="score",
                             high_line = 0.5,
                             group_var = "disease",
                             group_colors=c("healthy" = "steelblue",
                                            "CRC"= "brown3"),
                             point_size = 2,
                             bg_colors = NULL,
                             jitter_width = 0.1)
  expect_true(is.ggplot(p))
  expect_equal(colnames(p$data)[26], "sam.ind.x")

})
