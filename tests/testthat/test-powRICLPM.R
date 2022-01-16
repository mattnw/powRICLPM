test_that("à priori power analysis using powRICLPM() runs", {
  # Create valid powRICLPM() input
  Phi <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  wSigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2, byrow = TRUE)

  # Scenario 1: Inadmissible results and no recommended sample size
  test1 <- powRICLPM(target_power = 0.8,
                     parameter = "wB2~wA1",
                     search_lower = 50,
                     search_upper = 100,
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
  expect_equal(length(test1), 2)
  expect_equal(names(test1), c("conditions", "session"))
  expect_snapshot_output(powRICLPM_summary(test1, parameter = "wB2~wA1"))

  # Scenario 2: Scenario 1, but on 2 cores
  test2 <- powRICLPM(target_power = 0.8,
                     parameter = "wB2~wA1",
                     search_lower = 50,
                     search_upper = 100,
                     search_step = 50,
                     time_points = 3,
                     ICC = 0.5,
                     RI_cor = 0.3,
                     Phi = Phi,
                     wSigma = wSigma,
                     reps = 30,
                     cores = 2,
                     seed = 123456)

  # Run tests
  expect_equal(length(test2), 2)
  expect_equal(names(test2), c("conditions", "session"))
  expect_snapshot_output(powRICLPM_summary(test2, parameter = "wB2~wA1"))

  # Scenario 3: larger sample size with recommendation and multicore
  test3 <- powRICLPM(target_power = 0.8,
                     parameter = "wB2~wA1",
                     search_lower = 400,
                     search_upper = 600,
                     search_step = 50,
                     time_points = 4,
                     ICC = 0.5,
                     RI_cor = 0.3,
                     Phi = Phi,
                     wSigma = wSigma,
                     reps = 100,
                     cores = 5,
                     seed = 123456)

  # Run tests
  expect_equal(length(test3), 2)
  expect_equal(names(test3), c("conditions", "session"))
  expect_snapshot_output(powRICLPM_summary(test3, parameter = "wB2~wA1"))

  # Scenario 4: post hoc analysis for verification of previous results
 test4 <- powRICLPM(sample_size = c(400, 450, 500),
                    time_points = c(3, 4),
                    ICC = 0.5,
                    RI_cor = 0.3,
                    Phi = Phi,
                    wSigma = wSigma,
                    reps = 200,
                    cores = 6,
                    seed = 123456)

 # Run tests
 expect_equal(length(test4), 2)
 expect_equal(names(test4), c("conditions", "session"))
 expect_snapshot_output(powRICLPM_summary(test4, parameter = "wB2~wA1"))
})
