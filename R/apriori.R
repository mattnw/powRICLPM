apriori <- function(target_power,
                    parameter,
                    search_lower,
                    search_upper,
                    search_step,
                    time_points,
                    ICC,
                    RI_cor,
                    Phi,
                    wSigma,
                    Psi,
                    reps,
                    cores,
                    seed,
                    save_path) {
  # Get candidate sample sizes
  sample_size <- seq(search_lower, search_upper, search_step)

  # Create lavaan syntax for generating data
  pop_synt <- create_lavaan(time_points = time_points,
                            ICC = ICC,
                            RI_cor = RI_cor,
                            Phi = Phi,
                            wSigma = wSigma,
                            Psi = Psi,
                            syntax = T)

  # Create lavaan parameter table for population model
  pop_tab <- create_lavaan(time_points = time_points,
                           ICC = ICC,
                           RI_cor = RI_cor,
                           Phi = Phi,
                           wSigma = wSigma,
                           Psi = Psi)

  # Create lavaan syntax for estimating the model
  est_synt <- create_lavaan(time_points = time_points,
                            syntax = T)

  # Create lavaan parameter table for estimating model
  est_tab <- create_lavaan(time_points = time_points)

  # Create list with each element containing info for a single condition
  conditions <- purrr::pmap(list(sample_size,
                                 time_points,
                                 ICC,
                                 reps,
                                 pop_synt,
                                 list(pop_tab),
                                 est_synt,
                                 list(est_tab),
                                 save_path),
                            list)

  # Name all condition elements
  conditions <- purrr::map(conditions, function(x) {
    names(x) <- c("sample_size", "time_points", "ICC", "reps", "pop_synt", "pop_tab", "est_synt", "est_tab", "save_path")
    return(x)
  })

  # Create powRICLPM object, combining `conditions` and general session info
  object <- list(
    conditions = conditions,
    session = list(type = "apriori",
                   target_power = target_power,
                   parameter = parameter,
                   search_lower = search_lower,
                   search_upper = search_upper,
                   search_step = search_step,
                   reps = reps)
  )

  # Run Monte Carlo simulation for each condition
  if (cores == 1) {
    set.seed(seed)
    object$conditions <- purrr::map(object$conditions, run_condition)
  } else {
    future::plan(future::multisession, workers = cores)
    object$conditions <- furrr::future_map(object$conditions, run_condition,
                                .options = furrr::furrr_options(seed = seed))
    future::plan(future::sequential)
  }

  # Display results for a specific parameter
  summary.powRICLPM(object = object, parameter = parameter)

  # Make object of type `powRICLPM`
  class(object) <- "powRICLPM"
  return(object)
}


