test_that("Mplus model syntax creation works", {
  # Create valid powRICLPM_Mplus() input
  Phi <- matrix(c(.5, .1, .4, .5), ncol = 2, byrow = TRUE)
  wSigma <- matrix(c(1 , .3, .3, 1), ncol = 2, byrow = TRUE)
  Psi <- NULL
  save_path <- "./saved"

  # Create Mplus input file for RI-CLPM MCMC power analysis
  Mplus_table <- powRICLPM_Mplus(sample_size = 300,
                                 time_points = 4,
                                 ICC = 0.5,
                                 RI_cor = 0.3,
                                 Phi = Phi,
                                 wSigma = wSigma,
                                 reps = 10000,
                                 seed = 123456,
                                 save_path = save_path)

  # Test if directory and Mplus file exists
  expect_true(dir.exists(save_path), TRUE)
  expect_true(file.exists(file.path(save_path, paste0("Mplus_N", 300, "T", 4, ".txt"))))
})
