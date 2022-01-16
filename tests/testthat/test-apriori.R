test_that("apriori() works", {
  # Create valid apriori() input
  target_power <- 0.8
  parameter <- "wB2~wA1"
  search_lower <- 300
  search_upper <- 400
  search_step <- 50
  time_points <- 3
  ICC <- 0.5
  RI_cor <- 0.3
  Phi <- matrix(c(0.4, 0.15, 0.2, 0.3), ncol = 2, byrow = TRUE)
  wSigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2, byrow = TRUE)
  Psi <- compute_Psi(Phi = Phi, wSigma = wSigma)
  reps <- 30
  cores <- 1
  seed <- 123456
  save_path <- NA

  # Run apriori()
  output_apriori <- apriori(target_power = target_power,
                            parameter = parameter,
                            search_lower = search_lower,
                            search_upper = search_upper,
                            search_step = search_step,
                            time_points = time_points,
                            ICC = ICC,
                            RI_cor = RI_cor,
                            Phi = Phi,
                            wSigma = wSigma,
                            Psi = Psi,
                            reps = reps,
                            cores = cores,
                            seed = seed,
                            save_path = save_path)

  # Run general tests
  expect_equal(length(output_apriori), 2)
  expect_equal(names(output_apriori), c("conditions", "session"))

  # Test "conditions" element
  expect_equal(length(output_apriori$conditions), length(seq(search_lower, search_upper, search_step)))
  expect_type(output_apriori$conditions[[1]], "list")

  # Test "session" element
  expect_type(output_apriori$session, "list")
  expect_equal(length(output_apriori$session), 7)
  expect_equal(names(output_apriori$session), c("type", "target_power", "parameter", "search_lower", "search_upper", "search_step", "reps"))

})
