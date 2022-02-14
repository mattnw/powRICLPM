test_that("setup() works", {

  # Create valid setup() input
  target_power <- 0.8
  sample_size <- c(300, 350, 400)
  time_points <- 3
  ICC <- 0.5
  RI_cor <- 0.3
  Phi <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  wSigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2, byrow = TRUE)
  Psi <- compute_Psi(Phi = Phi, wSigma = wSigma)
  alpha <- 0.1
  reps <- 30
  seed <- 123456
  save_path <- NA

  # Run setup()
  output <- setup(target_power = target_power,
                  sample_size = sample_size,
                  time_points = time_points,
                  ICC = ICC,
                  RI_cor = RI_cor,
                  Phi = Phi,
                  wSigma = wSigma,
                  Psi = Psi,
                  skewness = 0,
                  kurtosis = 0,
                  alpha = alpha,
                  reps = reps,
                  seed = seed,
                  save_path = save_path)

  # Run general tests
  expect_equal(length(output), 2)
  expect_type(output, "list")
  expect_equal(names(output), c("conditions", "session"))

  # Test "conditions" element
  expect_equal(length(output$conditions), length(sample_size))
  expect_type(output$conditions[[1]], "list")

  # Test "session" element
  expect_type(output$session, "list")
  expect_equal(length(output$session), 6)
  expect_equal(names(output$session), c("sample_size", "time_points", "ICC", "reps", "target_power", "save_path"))

})
