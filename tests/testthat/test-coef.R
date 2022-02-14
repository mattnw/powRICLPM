test_that("coef_powRICLPM() works", {

  # Create valid powRICLPM() input
  Phi <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  wSigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2, byrow = TRUE)

  # Create powRICLPM object for à priori power analysis
  output <- powRICLPM(target_power = 0.5,
                      search_lower = 500,
                      search_upper = 600,
                      search_step = 50,
                      time_points = 3,
                      ICC = 0.5,
                      RI_cor = 0.3,
                      Phi = Phi,
                      wSigma = wSigma,
                      reps = 5,
                      seed = 123456)

  # Execute coef_powRICLPM()
  x <- coef_powRICLPM(output, parameter = "wB2~wA1")

  # Run tests
  expect_s3_class(x, "data.frame")
  expect_equal(dim(x), c(3, 14))
  expect_equal(names(x),
               c("sample_size", "time_points", "ICC", "errors",
                 "not_converged", "inadmissible", "pv", "avg",
                 "stdDev", "SEAvg", "mse", "acc", "cover", "pwr"))
})
