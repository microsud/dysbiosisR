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
  cloud.exp <- c(-0.02864206, -0.14630854, -0.08757229,
               -0.13069932, -0.18592467,  0.07201936)
  expect_equal(cloud.test,cloud.exp, tolerance=1e-3)

  # check score location
  expect_equal(colnames(cloud.results)[1:4], c("stats", "pvals", "ids","log2Stats"))
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
  cloud.exp2 <- c(-0.21524677, -0.35087536, -0.32974244,
                 -0.30326279, -0.47163967,  -0.07356939)
  expect_equal(cloud.test2,cloud.exp2, tolerance=1e-3)

  # check score location
  expect_equal(colnames(cloud.results2)[1:4], c("stats", "pvals", "ids","log2Stats"))
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
