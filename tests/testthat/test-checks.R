
# Test check_T() ----
test_that("check_T() works", {
  expect_equal(check_T(c(3, 4)), c(3, 4))
  expect_error(check_T(3.5), "Elements in `time_points` should be integers.")
  expect_error(check_T(c(2, 3)), "Elements in `time_points` should be larger than 2: The RI-CLPM is not identified with fewer than 3 time points.")
  expect_warning(check_T(c(3:30)), "A large number of repeated measures can lead to computational problems. You want to consider methods for intensive longitudinal data.")
})

# Test check_PD() ----
test_that("check_PD() works", {
  # Create sample matrices
  m1 <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
  m2 <- matrix(c(.3, .4, .3, .2), ncol = 2, byrow = TRUE)
  m3 <- matrix(c(1.2, .4, .3, 1.5), ncol = 2, byrow = TRUE)
  m4 <- matrix(c(.3, .4, -.3, .2), ncol = 2, byrow = TRUE)

  expect_equal(check_PD(m1, unit = TRUE), m1)
  expect_error(check_PD(m2, unit = TRUE), "Matrix is not positive definite.")
  expect_error(check_PD(m3, unit = TRUE), "Eigenvalues are not within unit circle.")
  expect_error(check_PD(m4, unit = TRUE), "Complex eigenvalues: Matrix is not positive definite.")
})

# Test check_N() ----
test_that("check_N() works", {
  expect_equal(check_N(c(200, 300), 3), c(200, 300))
  expect_error(check_N(c(200.4, 300), 3), "Elements in `sample_size` should be integers.")
  expect_error(check_N(c(-200, 300), 3), "sample_size` should be a positive integer.")
  expect_error(check_N(10, 3), "The number of parameters to be estimated is larger than the sample size.")
})

# Test check_RI() ----
test_that("check_ICC() works", {
  expect_equal(check_ICC(0.5), 0.5)
  expect_error(check_ICC(2), "`ICC` should be between 0 and 1.")
  expect_error(check_ICC("0.5"), "Elements in `ICC` should be of type `double`.")
})

# Test check_wSigma() ----
test_that("check_wSigma() works", {
  # Create sample matrices
  m1 <- matrix(c(1, .3, .3, 1), ncol = 2, byrow = TRUE)
  m2 <- matrix(c(.3, .4, .3, .2), ncol = 2, byrow = TRUE)
  m3 <- matrix(c(1, 2, 2, 1), ncol = 2, byrow = TRUE)

  expect_equal(check_wSigma(m1), m1)
  expect_error(check_wSigma(m2), "`wSigma` should be correlation matrix; 1's on the diagonal.")
  expect_error(check_wSigma("m2"), "`wSigma` should be of type `matrix`.")
  expect_error(check_wSigma(m3), "Matrix is not positive definite.")
})

# Test check_Phi() ----
test_that("check_Phi() works", {
  # Create sample matrices
  m1 <- matrix(c(.3, .2, .15, .2), ncol = 2, byrow = TRUE)
  m2 <- matrix(c(.8, .5, .4, .9), ncol = 2, byrow = TRUE)

  expect_equal(check_Phi(m1), m1)
  expect_error(check_Phi("m1"), "`Phi` should be of type `matrix`.")
  expect_error(check_Phi(m2), "Eigenvalues are not within unit circle.")
})

# Test check_skewness_kurtosis() ----
test_that("check_skewness_kurtosis() works", {
  expect_equal(check_skewness_kurtosis(c(0, 0.1, 0.6, 0.45)), c(0, 0.1, 0.6, 0.45))
  expect_error(check_skewness_kurtosis(c("0", 0.4)), "`skewness` and `kurtosis` should be of type `numeric`, `integer`, or `double`.")
})

# Test check_alpha() ----
test_that("check_alpha() works", {
  expect_equal(check_alpha(0.05), 0.05)
  expect_error(check_alpha("0.05"), "`alpha` should be of type `double`.")
  expect_error(check_alpha(-0.05), "`alpha` should be between 0 and 1.")
})

# Test check_seed() ----
test_that("check_seed() works", {
  expect_equal(check_seed(1234), 1234)
  expect_warning(check_seed(NA), "No seed was specified: A seed was randomly created.")
  expect_error(check_seed("1234"), "`seed` should be of type `numeric`.")
  expect_error(check_seed(1234.5), "`seed` should be an integer.")
})

# Test check_parameter() ----
test_that("check_parameter() works", {
  expect_equal(check_parameter("wB2~wA1"), "wB2~wA1")
  expect_error(check_parameter(1234), "`parameter` should be a character string specifying the parameter of interest.")
})

# Test check_reps() ----
test_that("check_reps() works", {
  expect_equal(check_reps(1000), 1000)
  expect_error(check_reps("1000"), "`reps` should be of type `numeric`.")
  expect_error(check_reps(1000.5), "`reps` should be a positive integer.")
  expect_error(check_reps(-1000), "`reps` should be a positive integer.")
})

# Test check_search() ----
test_that("check_search() works", {
  expect_error(check_search(100, 600, NULL), "`search_lower`, `search_upper`, and/or `search_step` were not declared. Please provided values for all `search_` arguments.")
  expect_invisible(check_search(100, 2000, 50))
  expect_error(check_search(2000, 100, 50), "`search_upper` should be higher than `search_lower`.")
  expect_error(check_search(100, 2000, 5000), "`search_step` should be smaller than or equal to the interval between `search_lower` and `search_upper`.")
  expect_error(check_search("100", 2000, 50), "`search_lower`, `search_upper`, and `search_step` should be of type `numeric`.")
})

# Test check_target() ----
test_that("check_target() works", {
  expect_equal(check_target(0.80), 0.80)
  expect_error(check_target("0.80"), "`target` should be of type `numeric`.")
  expect_error(check_target(1.4), "`target` should be between 0 and 1.")
})

# Test check_parameter_summary() ----
test_that("check_parameter_summary() works", {

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

  # Run test
  expect_snapshot_error(summary(output, parameter = "wB2~wB3"))
  })

