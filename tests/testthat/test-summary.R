test_that("summary.powRICLPM() works", {

  # Create powRICLPM object for à priori power analysis
  output <- powRICLPM(target_power = 0.5,
                      parameter = "wB2~wA1",
                      search_lower = 500,
                      search_upper = 600,
                      search_step = 50,
                      time_points = 3,
                      ICC = 0.5,
                      RI_cor = 0.3,
                      Phi = Phi,
                      wSigma = wSigma,
                      reps = 30,
                      cores = 1,
                      seed = 123456)

  # Run tests
  expect_error(summary(output), "Please provide a `parameter` argument. You can use the `names` function to get an overview of the parameter names in the powRICLPM object.")
  expect_error(summary(output, parameter = "wB2~wB3"), "The specified parameter is not valid. Use `names.powRICLPM` to get an overview of parameter names in the powRICLPM object.")
  expect_snapshot(cat(summary(output1, parameter = "wB2~wA1")))
})
