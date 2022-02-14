setup <- function(target_power,
                  sample_size,
                  time_points,
                  ICC,
                  RI_cor,
                  Phi,
                  wSigma,
                  Psi,
                  skewness,
                  kurtosis,
                  alpha,
                  reps,
                  seed,
                  save_path) {

  # Create sample_size, time_points, and ICC elements for powRICLPM condition list
  grid <- expand.grid(sample_size = sample_size,
                      time_points = time_points,
                      ICC = ICC)

  # Create lavaan syntax for generating data
  pop_synt <- purrr::pmap(list(grid$time_points, grid$ICC),
                          function(time_points, ICC) {
                            create_lavaan(time_points = time_points,
                                          ICC = ICC,
                                          RI_cor = RI_cor,
                                          Phi = Phi,
                                          wSigma = wSigma,
                                          Psi = Psi,
                                          syntax = TRUE)
                            }
                          )

  # Create lavaan parameter table for population model
  pop_tab <- purrr::pmap(list(grid$time_points, grid$ICC),
                         function(time_points, ICC) {
                           create_lavaan(time_points = time_points,
                                         ICC = ICC,
                                         RI_cor = RI_cor,
                                         Phi = Phi,
                                         wSigma = wSigma,
                                         Psi = Psi)
                         }
                         )

  # Create lavaan syntax for estimating the model
  est_synt <- purrr::map(grid$time_points, create_lavaan, syntax = TRUE)

  # Create lavaan parameter table for estimating model
  est_tab <- purrr::map(grid$time_points, create_lavaan)

  # Create list with each element containing info for a single condition
  conditions <- purrr::pmap(list(grid$sample_size,
                             grid$time_points,
                             grid$ICC,
                             reps,
                             pop_synt,
                             pop_tab,
                             est_synt,
                             est_tab,
                             skewness,
                             kurtosis,
                             alpha,
                             save_path),
                        list)

  # Name elements in powRICLPM condition list
  conditions <- purrr::map(conditions, function(x) {
    names(x) <- c("sample_size", "time_points", "ICC", "reps",
                  "pop_synt", "pop_tab", "est_synt", "est_tab",
                  "skewness", "kurtosis", "alpha",
                  "save_path")
    return(x)
  })

  # Create powRICLPM object, combining conditions and general session info
  object <- list(
    conditions = conditions,
    session = list(sample_size = sample_size,
                   time_points = time_points,
                   ICC = ICC,
                   reps = reps,
                   target_power = target_power,
                   save_path = save_path)
  )

  return(object)
}
