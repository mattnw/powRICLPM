posthoc <- function(sample_size,
                    time_points,
                    ICC,
                    RI_cor,
                    Phi,
                    Psi,
                    wSigma,
                    reps,
                    cores,
                    seed,
                    save_path) {

  # Create sample_size and time_points elements for powRICLPM condition list
  sample_size_list <- rep(sample_size, times = length(time_points))
  time_points <- rep(time_points, each = length(sample_size))

  # Create lavaan syntax for generating data
  pop_synt <- purrr::map(time_points, create_lavaan,
                         ICC = ICC,
                         RI_cor = RI_cor,
                         Phi = Phi,
                         wSigma = wSigma,
                         Psi = Psi,
                         syntax = TRUE)

  # Create lavaan parameter table for population model
  pop_tab <- purrr::map(time_points, create_lavaan,
                        ICC = ICC,
                        RI_cor = RI_cor,
                        Phi = Phi,
                        wSigma = wSigma,
                        Psi = Psi)

  # Create lavaan syntax for estimating the model
  est_synt <- purrr::map(time_points, create_lavaan, syntax = TRUE)

  # Create lavaan parameter table for estimating model
  est_tab <- purrr::map(time_points, create_lavaan)

  # Create list with each element containing info for a single condition
  conditions <- purrr::pmap(list(sample_size_list,
                             time_points,
                             ICC,
                             reps,
                             pop_synt,
                             pop_tab,
                             est_synt,
                             est_tab,
                             save_path),
                        list)

  # Name elements in powRICLPM condition list
  conditions <- purrr::map(conditions, function(x) {
    names(x) <- c("sample_size", "time_points", "ICC", "reps", "pop_synt", "pop_tab", "est_synt", "est_tab", "save_path")
    return(x)
  })

  # Create powRICLPM object, combining `conditions` and general session info
  object <- list(
    conditions = conditions,
    session = list(type = "posthoc",
                   sample_size = sample_size,
                   time_points = time_points,
                   ICC = ICC,
                   reps = reps)
  )

  # Run Monte Carlo simulation for each condition
  if(cores == 1) {
    set.seed(seed)
    object$conditions <- purrr::map(object$conditions, run_condition)
  } else {
    future::plan(future::multisession, workers = cores)
    object$conditions <- furrr::future_map(object$conditions, run_condition,
                                .options = furrr::furrr_options(seed = seed),
                                .progress = TRUE)
    future::plan(future::sequential)
  }

  # Display results for a specific parameter
  powRICLPM_summary(object = object)

  return(object)
}
