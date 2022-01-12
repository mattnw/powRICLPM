test_that("Mplus model syntax creation works", {
  # Create valid lav_...() input
  Phi <- matrix(c(.5, .1, .4, .5), ncol = 2, byrow = T)
  wSigma <- matrix(c(1 , .3, .3, 1), ncol = 2, byrow = T)
  Psi <- NULL

  # Create Mplus input file for RI-CLPM MCMC power analysis
  Mplus_table <- powRICLPM_Mplus(sample_size = 300,
                                 time_points = 4,
                                 ICC = 0.5,
                                 RI_cor = 0.3,
                                 Phi = Phi,
                                 wSigma = wSigma,
                                 reps = 10000,
                                 seed = 123456,
                                 save_path = "./saved")

  # Test if directory and Mplus file exists
  expect_true(dir.exists("./saved"), T)
  expect_true(file.exists(file.path(save_path, paste0("Mplus_N", sample_size, "T", time_points, ".txt"))))
})
