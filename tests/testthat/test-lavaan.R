test_that("lavaan model syntax creation works", {
  # Create valid create_lavaan() input
  Phi <- matrix(c(.5, .1, .4, .5), ncol = 2, byrow = T)
  wSigma <- matrix(c(1 , .3, .3, 1), ncol = 2, byrow = T)
  Psi <- matrix(c(0.71, -0.037, -0.037, 0.47), ncol = 2, byrow = T)

  # Generate lavaan parameter table
  pop_tab <- create_lavaan(time_points = 4,
                              ICC = 0.5,
                              RI_cor = 0.3,
                              Phi = Phi,
                              wSigma = wSigma,
                              Psi = Psi,
                              syntax = F)

  # Generate lavaan model syntax
  pop_synt <- create_lavaan(time_points = 3,
                              ICC = 0.5,
                              RI_cor = 0.3,
                              Phi = Phi,
                              wSigma = wSigma,
                              Psi = Psi,
                              syntax = T)

  # Create lavaan syntax for estimating a model
  est_synt <- create_lavaan(time_points = 3,
                            syntax = T)

  # Create lavaan parameter table for estimating model
  est_tab <- create_lavaan(time_points = 3)

  # Test parameter tables
  expect_type(pop_tab, "list")
  expect_type(est_tab, "list")

  expect_equal(dim(pop_tab), c(51, 5))
  expect_equal(dim(est_tab), c(38, 5))

  # Test syntax
  expect_type(pop_synt, "character")
  expect_equal(length(pop_synt), 1)

  expect_type(est_synt, "character")
  expect_equal(length(est_synt), 1)
})
