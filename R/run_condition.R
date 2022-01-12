#' @title
#' Run Monte Carlo Simulation For Single Condition
#'
#' @description
#' This function runs a Monte Carlo simulation for a given condition. It generates data based on the \code{pop_synt} element in "object", and estimates an RI-CLPM using the \code{est_synt`} element in "object". Data generation and model estimation are done using \pkg{lavaan}.
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
#'   \item \code{results}: A data frame containing the results (i.e., population values, bias, standard error of the estimate, 95% coverage, power, etc.), and
#'   \item \code{fail}: A logical vector denoting the Monte Carlo replications that resulted in negative variances, non-positive definite matrices, or fatal errors.
#' }
#'
#' @return A list.
#'
#' @importFrom lavaan simulateData inspect parameterEstimates coef
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
  sigs <- cover95 <- matrix(F, nrow = n_parameters, ncol = object$reps)

  # Initialize fail counter and index
  fails <- rep(F, times = object$reps)

  # Create folder for saving data
  if(!is.na(object$save_path)) {
    save_path <- file.path(
      object$save_path, paste0(
        "data_N", object$sample_size, "T", object$time_points, "ICC", object$ICC
      )
    )
    dir.create(save_path)
  }

  # Create lavaan wrapper that captures unintended side-effects (preventing fatal errors)
  possibly_lavaan <- purrr::possibly(lavaan::lavaan, otherwise = NA)

  # Start simulation
  for (r in 1:object$reps) {
    # Generate data
    dat <- simulateData(object$pop_synt, sample.nobs = object$sample_size)

    # Save data
    if(!is.na(object$save_path)) {
      fwrite(dat,
             file = file.path(object$save_path, paste0("df", r, ".dat")),
             sep = "\t", col.names = F, row.names = F, na = "-999")
    }

    # Fit model
    fit <- suppressWarnings(
      possibly_lavaan(object$est_synt, dat)
    )

    # Check negative variances, non-positive definite matrices, or fatal errors
    if (!inspect(fit, "post.check") || !isS4(fit)) {
      fails[r] <- T
      next
    }

    # Get estimates
    coefs[, r] <- coef(fit) # Save coefficients
    SEs[, r] <- parameterEstimates(fit, remove.nonfree = T)$se # Save standard errors
    sigs[, r] <- parameterEstimates(fit, remove.nonfree = T)$pvalue < .05
    low95[, r] <- parameterEstimates(fit, remove.nonfree = T)$ci.lower
    up95[, r] <- parameterEstimates(fit, remove.nonfree = T)$ci.upper
  }

  # Create and save repList
  if(!is.na(object$save_path)) {
    df_list <- paste0("df", 1:object$reps, ".dat")

    # Delete names of failed replications
    if(sum(fails) != 0) {
      df_list <- df_list[-which(fails == T)]
    }

    # Save list of generated data sets
    fwrite(as.list(df_list),
           file = file.path(save_path, "dfList.dat"),
           sep = "\n", col.names = F, row.names = F)
  }

  # Compute simulation results
  successes <- object$reps - sum(fails)
  avg <- rowMeans(coefs, na.rm = T)
  stdDev <- apply(coefs, 1, sd, na.rm = T)
  SEAvg <- rowMeans(SEs, na.rm = T)
  mse <- rowMeans((coefs - pv)^2, na.rm = T)
  cover95 <- rowSums(pv > low95 & pv < up95, na.rm = T) / successes
  sig <- rowSums(sigs, na.rm = T) / successes

  # Bind results of this particular scenario and add to object
  object$results <- data.frame(par, pv, avg, stdDev, SEAvg, mse, cover95, sig)
  object$fails <- fails

  return(object)
}
