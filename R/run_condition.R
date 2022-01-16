#' @title
#' Run Monte Carlo Simulation For Single Condition
#'
#' @description
#' This function runs a Monte Carlo simulation for a given condition. It generates data based on the \code{pop_synt} element in "object", and estimates an RI-CLPM using the \code{est_synt} element in "object". Data generation and model estimation are done using \pkg{lavaan}.
#'
#' @param object A list with information for running a single simulation condition. See "Details" for an overview of the elements this object must contain.
#'
#' @details
#' \subsection{Input: Elements in "object"}
#' To successfully run a Monte Carlo simulation, the "object" arguments needs the following elements:
#' \itemize{
#'   \item \code{sample_size}: The sample size.
#'   \item \code{time_points}: The number of time points.
#'   \item \code{reps}: Number of replications.
#'   \item \code{save_path}: The directory to save (data) files to.
#'   \item \code{pop_synt}: \pkg{lavaan} model syntax containing population values for data generation.
#'   \item \code{pop_tab}: \pkg{lavaan} parameter table for data generation.
#'   \item \code{est_synt}: \pkg{lavaan} model syntax for estimation.
#'   \item \code{est_tab}: \pkg{lavaan} parameter table for estimation.
#' }
#'
#' \subsection{Output}
#' This function adds the following two elements to "object":
#' \itemize{
#'   \item \code{results}: A data frame containing the results (i.e., population values, bias, standard error of the estimate, 95% coverage, power, etc.),
#'   \item \code{errors}: A logical vector denoting failed Monte Carlo replications,
#'   \item \code{not_converged}: A logical vector denoting non-converged Monte Carlo replications, and
#'   \item \code{inadmissible} A logical vector denoting Monte Carlo replications that resulted in negative variances or non-positive definite matrices.
#' }
#'
#' \subsection{Save}
#' This function has the option to save all generated datasets.
#'
#' @return A list.
#'
#' @importFrom lavaan simulateData inspect parameterEstimates coef lavaan
#' @importFrom data.table fwrite
run_condition <- function(object) {
  # Number of variables
  k <- 2

  # Compute number of parameters
  n_parameters <-
    sum(factorial(1 + k) / (2*(k - 1)) * (object$time_points + 1), # Number of (co)variances within- and between-level
        k^2*(object$time_points - 1)) # Number of lagged effects

  # Get index of freely estimated parameters
  index_parameter_free <- object$est_tab$pv == ""

  # Get the population values
  pv <- as.numeric(object$pop_tab$pv)[index_parameter_free]
  par <- paste0(object$pop_tab$lhs[index_parameter_free],
                object$pop_tab$op[index_parameter_free],
                object$pop_tab$rhs[index_parameter_free])

  # Allocate memory for results
  coefs <- SEs <- low95 <- up95 <- matrix(NA, nrow = n_parameters, ncol = object$reps)
  sigs <- cover95 <- matrix(FALSE, nrow = n_parameters, ncol = object$reps)

  # Initialize indeces for estimation problems
  errors <- warnings <- not_converged <- inadmissible <- rep(FALSE, times = object$reps)

  # Create folder for saving data
  if(!is.na(object$save_path)) {
    save_path <- file.path(
      object$save_path, paste0(
        "data_N", object$sample_size, "T", object$time_points, "ICC", object$ICC
      )
    )
    dir.create(save_path)
  }

  # Create lavaan function for safe and quiet error and warning handling
  safe_quiet_lavaan <- purrr::safely(purrr::quietly(lavaan::lavaan))
  quiet_lavInspect <- purrr::quietly(lavaan::lavInspect)

  # Start simulation
  for (r in 1:object$reps) {
    # Generate data
    dat <- simulateData(object$pop_synt, sample.nobs = object$sample_size)

    # Save data
    if(!is.na(object$save_path)) {
      fwrite(dat,
             file = file.path(object$save_path, paste0("df", r, ".dat")),
             sep = "\t", col.names = FALSE, row.names = FALSE, na = "-999")
    }

    # Fit model
    fit <- safe_quiet_lavaan(object$est_synt, data = dat)

    # Check if fatal error occurred
    if (!is.null(fit$error)) {
      errors[r] <- TRUE
      next
    }

    if (!identical(fit$result$warnings, character(0))) {

      # Check if solution converged
      if (!inspect(fit$result$result, what = "converged")) {
        not_converged[r] <- TRUE
      }

      # Check if parameter estimates are admissible
      if (quiet_lavInspect(fit$result$result, what = "post.check")$result != TRUE) {
        inadmissible[r] <- TRUE
      }
    }

    # Get estimates
    coefs[, r] <- coef(fit$result$result) # Save coefficients
    SEs[, r] <- parameterEstimates(fit$result$result, remove.nonfree = TRUE)$se # Save standard errors
    sigs[, r] <- parameterEstimates(fit$result$result, remove.nonfree = TRUE)$pvalue < .05
    low95[, r] <- parameterEstimates(fit$result$result, remove.nonfree = TRUE)$ci.lower
    up95[, r] <- parameterEstimates(fit$result$result, remove.nonfree = TRUE)$ci.upper
  }

  # Create and save repList
  if(!is.na(object$save_path)) {
    df_list <- paste0("df", 1:object$reps, ".dat")

    # Delete names of non-converged replications
    if(sum(not_converged) != 0) {
      df_list <- df_list[-which(not_converged)]
    }

    # Save list of generated data setss
    fwrite(as.list(df_list),
           file = file.path(save_path, "dfList.dat"),
           sep = "\n", col.names = FALSE, row.names = FALSE)
  }

  # Compute simulation results
  converged_reps <- object$reps - sum(not_converged)
  avg <- rowMeans(coefs, na.rm = TRUE)
  stdDev <- apply(coefs, 1, stats::sd, na.rm = TRUE)
  SEAvg <- rowMeans(SEs, na.rm = TRUE)
  mse <- rowMeans((coefs - pv)^2, na.rm = TRUE)
  cover95 <- rowSums(pv > low95 & pv < up95, na.rm = TRUE) / converged_reps
  sig <- rowSums(sigs, na.rm = TRUE) / converged_reps

  # Bind results of this particular scenario and add to object
  object$results <- data.frame(par, pv, avg, stdDev, SEAvg, mse, cover95, sig)
  object$errors <- errors
  object$not_converged <- not_converged
  object$inadmissible <- inadmissible

  return(object)
}





