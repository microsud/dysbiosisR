test_that("cloudStatistic works", {
  library(dysbiosisR)
  dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
  # get reference samples
  ref.samples <- sample_names(subset_samples(WirbelJ_2018,
                                             disease == "healthy"))
  cloud.results <- cloudStatistic(WirbelJ_2018,
                                  dist_mat = dist.mat,
                                  reference_samples = ref.samples,
                                  ndim=-1,
                                  k_num=80)
  cloud.test <- cloud.results$log2Stats[c(1:3,70:72)]
  cloud.exp <- c(-0.146593618, -0.046204044, -0.001462691,
                 0.025839734, -0.217112866, -0.199849374)
  expect_equal(cloud.test,cloud.exp, tolerance=1e-3)

  # check score location
  expect_equal(colnames(cloud.results)[1:3], c("stats", "pvals", "log2Stats"))
  #check all samples returned
  expect_equal(length(intersect(sample_names(WirbelJ_2018),
                                rownames(cloud.results))),
               length(sample_names(WirbelJ_2018)))
  # Change k num
  cloud.results2 <- cloudStatistic(WirbelJ_2018,
                                  dist_mat = dist.mat,
                                  reference_samples = ref.samples,
                                  ndim=-1,
                                  k_num=20)
  cloud.test2 <- cloud.results2$log2Stats[c(1:3,70:72)]
  cloud.exp2 <- c(-0.3547425, -0.3307632, -0.1813477,
                  -0.2242164, -0.4369376, -0.4498903)
  expect_equal(cloud.test2,cloud.exp2, tolerance=1e-3)

  # check score location
  expect_equal(colnames(cloud.results2)[1:3], c("stats", "pvals", "log2Stats"))
  #check all samples returned
  expect_equal(length(intersect(sample_names(WirbelJ_2018),
                                rownames(cloud.results2))),
               length(sample_names(WirbelJ_2018)))

  # check key inputs
  expect_error(cloudStatistic(WirbelJ_2018,
                                dist_mat = NULL,
                                reference_samples = ref.samples,
                                ndim=-1,
                                k_num=20))

  expect_error(cloudStatistic(WirbelJ_2018,
                              dist_mat = dist.mat,
                              reference_samples = NULL,
                              ndim=-1,
                              k_num=20))
})
