
# Test seed functionality with furrr package ----
test_that("setting a seed in furrr_options() works", {
  future::plan(future::multisession)
  x1 <- furrr::future_map(1:3, ~rnorm(10), .options = furrr::furrr_options(seed = 123))
  x2 <- furrr::future_map(1:3, ~rnorm(10), .options = furrr::furrr_options(seed = 123))
  future::plan(future::sequential)
  expect_equal(x1, x2)
}
)

# Test parameter order in lavaan package ----
test_that("lavaan() does not change order of parameters", {

  # Create valid powRICLPM() input
  time_points <- 3
  ICC <- 0.5
  RI_cor <- 0.3
  Phi <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  wSigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2, byrow = TRUE)
  Psi <- compute_Psi(Phi = Phi, wSigma = wSigma)

  # Create lavaan syntax for generating data
  pop_synt <- create_lavaan(time_points = time_points,
                            ICC = ICC,
                            RI_cor = RI_cor,
                            Phi = Phi,
                            wSigma = wSigma,
                            Psi = Psi,
                            syntax = TRUE)

  # Create lavaan parameter table for population model
  pop_tab <- create_lavaan(time_points = time_points,
                           ICC = ICC,
                           RI_cor = RI_cor,
                           Phi = Phi,
                           wSigma = wSigma,
                           Psi = Psi)

  # Create lavaan syntax for estimating the model
  est_synt <- create_lavaan(time_points = time_points,
                            syntax = TRUE)

  # Create lavaan syntax for estimating the model
  est_tab <- create_lavaan(time_points = time_points)

  # Get index free parameters
  index_parameter_free <- est_tab$pv == ""

  # Create parameter names from population model
  pop_par <- paste0(pop_tab$lhs[index_parameter_free],
                    pop_tab$op[index_parameter_free],
                    pop_tab$rhs[index_parameter_free])

  # Generate data
  dat <- simulateData(pop_synt)

  # Fit model
  fit <- lavaan(est_synt, data = dat)

  # Create parameter names from fitted model
  fit_par <- paste0(lavaan::parameterestimates(fit)$lhs[index_parameter_free],
                    lavaan::parameterestimates(fit)$op[index_parameter_free],
                    lavaan::parameterestimates(fit)$rhs[index_parameter_free])

  # Check if parameters are the same order
  expect_equal(pop_par, fit_par)
  }
)
