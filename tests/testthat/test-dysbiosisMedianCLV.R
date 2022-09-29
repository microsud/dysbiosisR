test_that("dysbiosisMedianCLV works", {
  library(dysbiosisR)
  dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
  # get reference samples
  ref.samples <- sample_names(subset_samples(WirbelJ_2018,
                                             disease == "healthy"))

  dysbiosis_1 <- dysbiosisMedianCLV(WirbelJ_2018,
                                    dist_mat = dist.mat,
                                    reference_samples = ref.samples)

  dys.test <- dysbiosis_1$score[c(1:3,70:72)]
  dys.exp <- c(0.7698956,
               0.7611880,
               0.6879019,
               0.5607490,
               0.6468704,
               0.6725262)

  expect_equal(dys.test,dys.exp, tolerance=1e-3)

  #check all samples returned
  expect_equal(length(intersect(sample_names(WirbelJ_2018),
                                rownames(dysbiosis_1))),
               length(sample_names(WirbelJ_2018)))

  # check score location
  expect_equal(colnames(dysbiosis_1)[1:3], c("score",
                                             "study_name",
                                             "subject_id"))
  # Check key inputs
  expect_error(dysbiosisMedianCLV(WirbelJ_2018,
                                  dist_mat = NULL,
                                  reference_samples = ref.samples))

  expect_error(dysbiosisMedianCLV(WirbelJ_2018,
                                  dist_mat = dist.mat,
                                  reference_samples = NULL))
})
