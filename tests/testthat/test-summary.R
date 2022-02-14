test_that("powRICLPM_summary() works", {

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

  # Run tests
  expect_snapshot(cat(summary(output, parameter = "wB2~wA1")))
})
