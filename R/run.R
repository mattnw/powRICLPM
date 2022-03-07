#' @title
#' Run Monte Carlo Simulation For Single Condition
#'
#' @description
#' \code{run_condition()} runs a Monte Carlo simulation for a single given condition. It generates data based on the \code{pop_synt} element in "object", and estimates an RI-CLPM using the \code{est_synt} element in "object". Data generation and model estimation are done using \pkg{lavaan}.
#'
#' @param object A list with information for running a single simulation condition. See "Details" for an overview of the elements this object must contain.
#' @param p A \pkg{progressr} object
#'
#' @details
#' \subsection{Input: Elements in "object"}{To successfully run a Monte Carlo simulation, the "object" arguments needs the following elements:
#' \itemize{
#'   \item \code{sample_size}: The sample size.
#'   \item \code{time_points}: The number of time points.
#'   \item \code{reps}: Number of replications.
#'   \item \code{pop_synt}: \pkg{lavaan} model syntax containing population values for data generation.
#'   \item \code{pop_tab}: \pkg{lavaan} parameter table for data generation.
#'   \item \code{est_synt}: \pkg{lavaan} model syntax for estimation.
#'   \item \code{est_tab}: \pkg{lavaan} parameter table for estimation.
#'   \item \code{skewness}: The skewness value(s) for the observed variables.
#'   \item \code{kurtosis}: The kurtosis value(s) for the observed variables.
#'   \item \code{alpha}: The significance criterion.
#'   \item \code{save_path}: Folder to which simulated data are saved (optional).
#'   }
#' }
#'
#' \subsection{Output}{This function adds the following elements to "object":
#' \itemize{
#'   \item \code{results}: A data frame containing the results (i.e., population values, bias, standard error of the estimate, coverage, power, etc.),
#'   \item \code{errors}: A logical vector denoting failed Monte Carlo replications,
#'   \item \code{not_converged}: A logical vector denoting non-converged Monte Carlo replications, and
#'   \item \code{inadmissible} A logical vector denoting Monte Carlo replications that resulted in negative variances or non-positive definite matrices.
#'   }
#' }
#'
#' @return A list.
#'
#' @importFrom lavaan simulateData inspect parameterEstimates coef lavaan
run_condition <- function(object, p) {

  # Number of variables
  k <- 2 # Fixed in v0.2.0-alpha and earlier

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
  coefs <- SEs <- low <- up <- matrix(NA, nrow = n_parameters, ncol = object$reps)
  sigs <- cover <- matrix(FALSE, nrow = n_parameters, ncol = object$reps)

  # Initialize index for estimation problems
  errors <- warnings <- not_converged <- inadmissible <- rep(FALSE, times = object$reps)

  # Create lavaan function for safe and quiet error and warning handling
  safe_quiet_lavaan <- purrr::safely(purrr::quietly(lavaan::lavaan))
  quiet_lavInspect <- purrr::quietly(lavaan::lavInspect)

  # Create folder for saving simulated data to (optional)
  if (!is.na(object$save_path)) {

    # Create path to condition-specific folder
    save_path_aux <- file.path(
      object$save_path,
      paste0("data_N", object$sample_size, "_T", object$time_points, "_ICC", object$ICC)
    )

    # Create condition folder
    dir.create(save_path_aux)
  }

  # Start simulation
  for (r in 1:object$reps) {

    # Generate data
    dat <- simulateData(model = object$pop_synt,
                        sample.nobs = object$sample_size,
                        skewness = object$skewness,
                        kurtosis = object$kurtosis)

    # Save data
    if (!is.na(object$save_path)) {
      utils::write.table(dat,
                         file = file.path(save_path_aux, paste0("df", r, ".dat")),
                         sep = "\t", col.names = FALSE, row.names = FALSE, na = "-999")
    }

    # Fit model
    fit <- safe_quiet_lavaan(object$est_synt, data = dat, estimator = "DWLS")

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
    coefs[, r] <- coef(fit$result$result)
    SEs[, r] <- parameterEstimates(fit$result$result, remove.nonfree = TRUE)$se
    sigs[, r] <- parameterEstimates(fit$result$result, remove.nonfree = TRUE)$pvalue < object$alpha
    low[, r] <- parameterEstimates(fit$result$result,
                                     remove.nonfree = TRUE,
                                     level = object$alpha)$ci.lower
    up[, r] <- parameterEstimates(fit$result$result,
                                    remove.nonfree = TRUE,
                                    level = object$alpha)$ci.upper
  }

  # Create and save repList
  if (!is.na(object$save_path)) {
    df_list <- paste0("df", 1:object$reps, ".dat")

    # Delete names of non-converged replications
    if(sum(not_converged) != 0) {
      df_list <- df_list[-which(not_converged)]
    }

    # Save list of generated data setss
    utils::write.table(as.list(df_list),
                       file = file.path(save_path_aux, "dfList.dat"),
                       sep = "\n", col.names = FALSE, row.names = FALSE)
  }

  # Compute simulation results
  converged_reps <- object$reps - sum(not_converged)
  avg <- rowMeans(coefs, na.rm = TRUE)
  stdDev <- apply(coefs, 1, stats::sd, na.rm = TRUE)
  SEAvg <- rowMeans(SEs, na.rm = TRUE)
  mse <- rowMeans((coefs - pv)^2, na.rm = TRUE)
  acc <- rowSums(up - low, na.rm = TRUE) / converged_reps
  cover <- rowSums(pv > low & pv < up, na.rm = TRUE) / converged_reps
  pwr <- rowSums(sigs, na.rm = TRUE) / converged_reps

  # Bind results of this particular scenario and add to object
  object$results <- data.frame(par, pv, avg, stdDev, SEAvg, mse, acc, cover, pwr)
  object$errors <- errors
  object$not_converged <- not_converged
  object$inadmissible <- inadmissible

  # Signal progress bar
  p()

  return(object)
}





