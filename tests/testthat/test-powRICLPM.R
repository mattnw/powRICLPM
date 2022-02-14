test_that("sample size recommendation using powRICLPM() works", {

  # Create valid powRICLPM() input
  Phi <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  wSigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2, byrow = TRUE)

  # Scenario 1: Inadmissible results and no recommended sample size
  test1 <- powRICLPM(target_power = 0.8,
                     search_lower = 50,
                     search_upper = 100,
                     search_step = 50,
                     time_points = 3,
                     ICC = 0.5,
                     RI_cor = 0.3,
                     Phi = Phi,
                     wSigma = wSigma,
                     reps = 5,
                     seed = 123456,
                     save_path = "/Users/jeroenmulder/SURFdrive/R packages/powRICLPM/tests/powRICLPM")

  # Run tests
  expect_equal(names(test1), c("conditions", "session"))
  expect_equal(class(test1), c("powRICLPM", "list"))
  expect_snapshot_output(summary(test1, parameter = "wB2~wA1"))

  # Scenario 2: Scenario 1, but in multisession using `furrr`
  future::plan(future::multisession, workers = 3)
  test2 <- powRICLPM(target_power = 0.8,
                     search_lower = 50,
                     search_upper = 100,
                     search_step = 50,
                     time_points = 3,
                     ICC = 0.5,
                     RI_cor = 0.3,
                     Phi = Phi,
                     wSigma = wSigma,
                     reps = 5,
                     seed = 123456)
  future::plan(future::sequential)

  # Run tests
  expect_equal(names(test2), c("conditions", "session"))
  expect_equal(class(test1), c("powRICLPM", "list"))
  expect_snapshot_output(summary(test2, parameter = "wB2~wA1"))

  # Scenario 3: Get performance measures
  future::plan(future::multisession, workers = 3)
  test3 <- powRICLPM(target_power = 0.8,
                     sample_size = c(400, 500),
                     time_points = c(3, 4),
                     ICC = c(0.3, 0.7),
                     RI_cor = 0.3,
                     Phi = Phi,
                     wSigma = wSigma,
                     reps = 5,
                     seed = 123456)
  future::plan(future::sequential)

 # Run tests
 expect_equal(names(test3), c("conditions", "session"))
 expect_equal(class(test1), c("powRICLPM", "list"))
 expect_snapshot_output(summary(test3, parameter = "wB2~wA1"))
})
