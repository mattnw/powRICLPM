#' @method summary powRICLPM
#' @export
summary.powRICLPM <- function(object, parameter = NULL, ...) {

  # Argument validation
  check_object(object)
  if (!is.null(parameter)) {parameter <- check_parameter_summary(parameter, object)}

  # Get target power
  target_power <- object$session$target_power

  # Print basic summary()-output
  cat("\n", "powRICLPM analysis completed:", sep = "")
  cat("\n", "- Monte Carlo replications: ", object$session$reps)
  cat("\n", "- Sample size(s):", object$session$sample_size)
  cat("\n", "- Number of time points:", object$session$time_points)
  cat("\n", "- Proportion(s) random intercept variance:", object$session$ICC)
  cat("\n", " - Target power: ", target_power, sep = "")

  if (!is.null(parameter) && !is.null(target_power)) {

    # Compute no conditions that reach targeted power
    n_recommendations <- length(
      purrr::keep(object$conditions, function(condition) {
        condition$results[condition$results$par == parameter, "pwr"] >= target_power
      })
    )

    # Print n_recommendations
    cat("\n\n", "Number of conditions that reach target power: ", n_recommendations, sep = "")

    if (n_recommendations == 0) {
    cat("\n\n", "Suggested next steps:", sep = "")
    cat("\n", "- Increase search_upper and rerun the analysis.")

    } else if (n_recommendations == 1) {

      # Detect condition that meets target_power
      candidate <- purrr::detect_index(object$conditions, function(condition) {
        condition$results[condition$results$par == parameter, "pwr"] >= target_power
      })

      # Print recommended experimental condition
      cat("\n", "- Sample size:", object$conditions[[candidate]]$sample_size)
      cat("\n", "- Number of time points:", object$conditions[[candidate]]$time_points)
      cat("\n", "- Proportion of between-unit variance:", object$conditions[[candidate]]$ICC)
      cat("\n", "- Number of non-converged replications:", sum(object$conditions[[candidate]]$not_converged))
      cat("\n", " - Number of replications with inadmissible estimates:", sum(object$conditions[[candidate]]$inadmissible), sep = "")

      cat("\n\n", "Sugested next step:", sep = "")
      cat("\n", "- If this is a preliminary powRICLPM analysis, validate this result by rerunning the analysis with an increased number of replications (e.g., `reps = 1000`).")

    } else if (n_recommendations > 1) {

      # Print suggested next steps
      cat("\n\n", "Suggested next steps:", sep = "")
      cat("\n", "- If this is a preliminary powRICLPM analysis, validate these recommendations (or a selection) by rerunning the analysis with an increased number of replications (e.g., `reps = 1000`).")
    }

    # Print suggestion that applies to every scenario
    cat("\n", "- Use `plot_powRICLPM()` to visualize results across all experimental conditions.")
  } else {

    # Print suggestions when no `parameter` was specified
    cat("\n\n", "Sugested next steps:", sep = "")
    cat("\n", " - Specify the `parameter` argument of `summary()` to obtain a parameter-specific summary.", sep = "")
    cat("\n", " - Use `names_powRICLPM()` to obtain parameter names in the powRICLPM object.", sep = "")
    cat("\n", " - Use `plot_powRICLPM()` to visualize results for a specific parameter across conditions.")
  }
}

#' @title
#' Performance Measures From `powRICLPM` Object
#'
#' @description
#' \code{coef_powRICLPM} extracts performance measures (e.g., bias, mean square error, power) for a specific parameter, across all experimental conditions, from a `powRICLPM` object.
#'
#' @param object A `powRICLPM` object.
#' @param parameter A character string denoting a single variable of interest. Use \code{\link{names_powRICLPM}} to get an overview of parameter names in the `powRICLPM` object.
#'
#' @return
#' A `data.frame` object with columns containing performance measures, and rows representing experimental conditions.
#' @export
coef_powRICLPM <- function(object, parameter) {

  # Argument validation
  check_object(object)
  parameter <- check_parameter_summary(parameter, object)

  # Combine sample sizes and simulated power across conditions
  d <- purrr::map_dfr(object$conditions, function(condition) {

    # Create data frame
    data.frame(sample_size = condition$sample_size,
               time_points = condition$time_points,
               ICC = condition$ICC,
               errors = sum(condition$errors),
               not_converged = sum(condition$not_converged),
               inadmissible = sum(condition$inadmissible),
               condition$results[condition$results$par == parameter, -1]
    )
  })
  return(d)
}


#' @title
#' Parameter Names From `powRICLPM` Object
#'
#' @description
#' \code{names_powRICLPM} extracts the names of the variables that are internally created by the powRICLPM package. Details about the naming conventions can be found in the "Details" section of \code{\link{powRICLPM}}.
#'
#' @param object A `powRICLPM` object.
#'
#' @details
#' When simulating the power for conditions with a varying number of time points, there are different amounts of parameters across the conditions. By default, this function returns the parameter names from the condition with the smallest number of parameters, such that the returned parameter names are valid for each condition.
#'
#' @return A character vector with the names of the variables internally created by the \pkg{powRICLPM} package.
#'
#' @examples
#' # Define population parameters for lagged effects and within-component correlations
#' Phi <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#' wSigma <- matrix(c(1, .3, .3, 1), ncol = 2, byrow = TRUE)
#'
#' # Create powRICLPM object for à priori power analysis
#' output <- powRICLPM(target_power = 0.5,
#'                     search_lower = 300,
#'                     search_upper = 600,
#'                     search_step = 50,
#'                     time_points = 3,
#'                     ICC = 0.5,
#'                     RI_cor = 0.3,
#'                     Phi = Phi,
#'                     wSigma = wSigma,
#'                     reps = 50,
#'                     seed = 123456)
#'
#' # Get names of internally created parameters
#' names_powRICLPM(output)
#' @export
names_powRICLPM <- function(object) {

  # Determine number of parameters per condition
  condition_length <- purrr::map_int(object$conditions, function(condition) {
    length(condition$results$par)
  })

  # Return names of parameters in condition with least parameters
  return(object$conditions[[which.min(condition_length)]]$results$par)
}


#' @method print powRICLPM
#' @export
print.powRICLPM <- function(x, ...) {
  cat("\n", "A powRICLPM object resulting from a call to powRICLPM():" , sep = "")
  cat("\n", "- Monte Carlo replications: ", x$session$reps)
  cat("\n", "- Sample size(s):", x$session$sample_size)
  cat("\n", "- Number of time points:", x$session$time_points)
  cat("\n", "- Proportion(s) random intercept variance:", x$session$ICC)
  cat("\n", "- Target power: ", x$sessions$target_power)
  invisible(x)
}
