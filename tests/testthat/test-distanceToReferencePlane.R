test_that("distanceToReferencePlane works", {
  library(dysbiosisR)
  dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
  # get reference samples
  ref.samples <- sample_names(subset_samples(WirbelJ_2018,
                                             disease == "healthy"))

  dtrp.results <- distanceToReferencePlane(WirbelJ_2018,
                                           dist_mat = dist.mat,
                                           reference_samples = ref.samples)

  dtrp.test <- dtrp.results$dtrpScore[c(3:6, 80:82)]
  dtrp.exp <- c(0.14400941,
                0.02996153,
                0.03891631,
                0.13630248,
                0.24036553,
                0.04886236,
                0.18372597)

  expect_equal(dtrp.test, dtrp.exp, tolerance = 1e-3)

  #check all samples returned
  expect_equal(length(intersect(sample_names(WirbelJ_2018),
                                rownames(dtrp.results))),
               length(sample_names(WirbelJ_2018)))

  # check score location
  expect_equal(colnames(dtrp.results)[1:3], c("dtrpScore",
                                             "study_name",
                                             "subject_id"))
  # Check key inputs
  expect_error(distanceToReferencePlane(WirbelJ_2018,
                                        dist_mat = NULL,
                                        reference_samples = ref.samples))

  expect_error(distanceToReferencePlane(WirbelJ_2018,
                                        dist_mat = dist.mat,
                                        reference_samples = NULL))
})
