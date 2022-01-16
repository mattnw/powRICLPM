test_that("run_condition() works", {
  # Create valid create_lavaan() input
  time_points <- 3
  ICC <- 0.5
  RI_cor <- 0.3
  Phi <- matrix(c(.5, .1, .4, .5), ncol = 2, byrow = TRUE)
  wSigma <- matrix(c(1 , .3, .3, 1), ncol = 2, byrow = TRUE)
  Psi <- matrix(c(0.710, -0.037, -0.037, 0.470), ncol = 2, byrow = TRUE)

  # Create lavaan syntax for generating data
  pop_synt <- create_lavaan(time_points = time_points,
                            ICC = ICC,
                            Phi = Phi,
                            Psi = Psi,
                            wSigma = wSigma,
                            RI_cor = RI_cor,
                            syntax = TRUE)

  # Create lavaan parameter table for population model
  pop_tab <- create_lavaan(time_points = time_points,
                           ICC = ICC,
                           Phi = Phi,
                           Psi = Psi,
                           wSigma = wSigma,
                           RI_cor = RI_cor)

  # Create lavaan syntax for estimating the model
  est_synt <- create_lavaan(time_points = time_points,
                            syntax = TRUE)

  # Create lavaan parameter table for estimating model
  est_tab <- create_lavaan(time_points = time_points)

  # Set condition information
  input <- list(sample_size = 300,
                time_points = 3,
                reps = 100,
                save_path = NA,
                pop_synt = pop_synt,
                pop_tab = pop_tab,
                est_synt = est_synt,
                est_tab = est_tab)

  # Run run_condition()
  output <- run_condition(object = input)

  # Test output
  expect_equal(length(output), length(input) + 4)
  expect_equal(names(output), c("sample_size", "time_points", "reps", "save_path", "pop_synt", "pop_tab", "est_synt", "est_tab", "results", "errors", "not_converged", "inadmissible"))
  expect_equal(length(output$inadmissible), input$reps)
  expect_type(output$results, "list")
  expect_type(output$errors, "logical")
  expect_type(output$not_converged, "logical")
  expect_type(output$inadmissible, "logical")
  expect_equal(dim(output$results), c(20, 8))
})
