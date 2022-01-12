test_that("compute_Psi() works", {
  # Set lagged effects
  Phi = matrix(c(.2, .15, .10, .3), ncol = 2, byrow = T)

  # Set variance-covariance of within-components
  wSigma = matrix(c(1 , .3, .3, 1) , ncol = 2, byrow = T)

  # Compute residual (co)variances
  output <- compute_Psi(Phi = Phi, wSigma = wSigma)

  # Run tests
  expect_type(output, "double")
  expect_equal(dim(output), c(2, 2))
  expect_equal(eigen(output)$values > 0, c(T, T))
})
