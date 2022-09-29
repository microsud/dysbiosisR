test_that("multiplication works", {
  library(dysbiosisR)
  # get reference samples
  ref.samples <- sample_names(subset_samples(WirbelJ_2018,
                                             disease == "healthy"))
  dysbiosis_3 <- combinedShannonJSD(WirbelJ_2018,
                                    reference_samples = ref.samples)

  sdjsd.test <- dysbiosis_3$ShannonJSDScore[c(1:3,70:72)]
  sdjsd.exp <- c(0.1988758761, 0.0768259552,  0.1964672414,
                 0.2118476444, -0.0004954855,  0.1109020483)
  expect_equal(sdjsd.test,sdjsd.exp, tolerance=1e-3)

  # check score location
  expect_equal(colnames(dysbiosis_3)[1], "ShannonJSDScore")
  #check all samples returned
  expect_equal(length(intersect(sample_names(WirbelJ_2018),
                                rownames(dysbiosis_3))),
               length(sample_names(WirbelJ_2018)))

  expect_equal(length(intersect(sample_names(WirbelJ_2018),
                                rownames(dysbiosis_3))),
               length(sample_names(WirbelJ_2018)))

  # check key inputs
  expect_error(combinedShannonJSD(WirbelJ_2018,
                                  reference_samples = NULL))

})
