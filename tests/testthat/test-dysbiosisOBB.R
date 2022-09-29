test_that("dysbiosisOBB", {
  library(dysbiosisR)
  dysbiosis.oob <- dysbiosisOBB(WirbelJ_2018,
                                group_col = "disease",
                                case_label = "CRC",
                                seed_value = 1235,
                                add_tuneRF_params = list(ntreeTry=100,
                                                         stepFactor=1.5,
                                                         improve=0.01,
                                                         trace=TRUE,
                                                         dobest=FALSE),
                                ntree = 100,
                                plot_roc = FALSE)
  obb.test <- dysbiosis.oob$oob.score[c(1:3,60:63)]
  obb.exp <- c(0.9459459, 0.7714286, 0.6969697,
               0.1578947, 0.2142857, 0.5405405, 0.2500000)
  expect_equal(obb.test,obb.exp, tolerance=1e-3)

  # check score location
  expect_equal(colnames(dysbiosis.oob)[1], "oob.score")
  #check all samples returned
  expect_equal(length(intersect(sample_names(WirbelJ_2018),
                               rownames(dysbiosis.oob))),
              length(sample_names(WirbelJ_2018)))

  expect_error(dysbiosisOBB(WirbelJ_2018,
                            group_col = NULL,
                            case_label = "CRC",
                            seed_value = 1235,
                            add_tuneRF_params = list(ntreeTry=100,
                                                     stepFactor=1.5,
                                                     improve=0.01,
                                                     trace=TRUE,
                                                     dobest=FALSE),
                            ntree = 100,
                            plot_roc = FALSE))

  expect_error(dysbiosisOBB(WirbelJ_2018,
                            group_col = "disease",
                            case_label = NULL,
                            seed_value = 1235,
                            add_tuneRF_params = list(ntreeTry=100,
                                                     stepFactor=1.5,
                                                     improve=0.01,
                                                     trace=TRUE,
                                                     dobest=FALSE),
                            ntree = 100,
                            plot_roc = FALSE))

})
