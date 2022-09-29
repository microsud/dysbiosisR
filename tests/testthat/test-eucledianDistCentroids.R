test_that("euclideanDistCentroids works", {

  library(dysbiosisR)
  dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
  dysbiosis_2 <- euclideanDistCentroids(WirbelJ_2018,
                                        dist_mat = dist.mat,
                                        use_squared = TRUE,
                                        group_col = "disease",
                                        control_label = "healthy",
                                        case_label = "CRC")

  cent.test <- dysbiosis_2$CentroidDist_score[c(1:3,70:72)]
  cent.exp <- c(0.07722225, 0.06710710, 0.05249821,
                -0.02977273, -0.01054721, -0.03099986)
  expect_equal(cent.test,cent.exp, tolerance=1e-3)

  # check score location
  expect_equal(colnames(dysbiosis_2)[1:3], c("CentroidDist_CRC",
                                             "CentroidDist_healthy",
                                             "CentroidDist_score"))
  #check all samples returned
  expect_equal(length(intersect(sample_names(WirbelJ_2018),
                                rownames(dysbiosis_2))),
               length(sample_names(WirbelJ_2018)))

  expect_error(euclideanDistCentroids(WirbelJ_2018,
                                      dist_mat = NULL,
                                      use_squared = TRUE,
                                      group_col = "disease",
                                      control_label = "healthy",
                                      case_label = "CRC"))

  expect_error(euclideanDistCentroids(WirbelJ_2018,
                                      dist_mat = dist.mat,
                                      use_squared = TRUE,
                                      group_col = NULL,
                                      control_label = "healthy",
                                      case_label = "CRC"))

  expect_error(euclideanDistCentroids(WirbelJ_2018,
                                      dist_mat = dist.mat,
                                      use_squared = TRUE,
                                      group_col = "disease",
                                      control_label = NULL,
                                      case_label = "CRC"))
  expect_error(euclideanDistCentroids(WirbelJ_2018,
                                      dist_mat = dist.mat,
                                      use_squared = TRUE,
                                      group_col = "disease",
                                      control_label = "healthy",
                                      case_label = NULL))

})

