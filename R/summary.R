#' @title
#' powRICLPM Power Analysis Summary
#'
#' @description
#' \code{powRICLPM_summary()} is a method for class powRICLPM, summarizing the results of an à priori or post hoc powRICLPM simulation study.
#'
#' @param object A powRICLPM object.
#' @param ... Phantom argument that has no influence on the summary produced.
#' @param parameter A character string denoting a single variable of interest. Use \code{\link{powRICLPM_names}} to get an overview of valid parameter names for the powRICLPM object.
#'
#' @seealso
#' \code{\link{powRICLPM_names}} to get an overview of valid parameter names included in the powRICLPM object.
#'
#' @examples
#' # Define population parameters for lagged effects and within-component correlations
#' Phi <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#' wSigma <- matrix(c(1, .3, .3, 1), ncol = 2, byrow = TRUE)
#'
#' # Create powRICLPM object for à priori power analysis
#' output <- powRICLPM(target_power = 0.5,
#'                      parameter = "wB2~wA1",
#'                      search_lower = 300,
#'                      search_upper = 600,
#'                      search_step = 50,
#'                      time_points = 3,
#'                      ICC = 0.5,
#'                      RI_cor = 0.3,
#'                      Phi = Phi,
#'                      wSigma = wSigma,
#'                      reps = 50,
#'                      cores = 1,
#'                      seed = 123456)
#'
#' # Summarize results
#' powRICLPM_summary(output, parameter = "wB2~wA1")
#' @export
powRICLPM_summary <- function(object, ..., parameter = NULL) {

  if (object$session$type == "apriori") {

    # Check arguments
    parameter = check_parameter_summary(parameter, object)

    # Get analysis characteristics
    target_power <- object$session$target_power

    # Find index recommended sample size
    candidate <- purrr::detect_index(object$conditions, function(x) {
      x$results[x$results$par == parameter, "sig"] >= target_power
    })

    # Display results
    cat("\n", "A priori powRICLPM-analysis completed:", sep = "")
    cat("\n", "- Monte Carlo replications: ", object$session$reps)
    cat("\n", "- Target power: ", target_power)

    if(candidate == 0) { # No recommended sample size found
      cat("\n", "- Preliminary recommended sample size: NOT FOUND")

      cat("\n\n", "Sugested next steps:", sep = "")
      cat("\n", "- Increase search_upper and rerun the analysis.")
      cat("\n", "- Use powRICLPM_plot() to visualize preliminary results for entire range of sample sizes.")

    } else { # Recommended sample size found

      cat("\n", "- Preliminary recommended sample size: ", object$conditions[[candidate]]$sample_size)
      cat("\n", "- Number of non-converged replications for recommended run: ", sum(object$conditions[[candidate]]$not_converged))
      cat("\n", "- Number of replications with inadmissible estimates for recommended run: ", sum(object$conditions[[candidate]]$inadmissible))

      cat("\n\n", "Sugested next steps:", sep = "")
      cat("\n", "- Use powRICLPM_plot() to visualize preliminary results for entire range of sample sizes.")
      cat("\n", "- Rerun powRICLPM with the recommended sample size and large number of replications to validate preliminary results.")
    }

  } else if (object$session$type == "posthoc") {

    if (is.null(parameter)) {

      # Display general results
      cat("\n", "Post hoc powRICLPM-analysis completed:", sep = "")
      cat("\n", "- Monte Carlo replications: ", object$session$reps)
      cat("\n", "- Proportion random intercept variance:", object$session$ICC)
      cat("\n", "- Sample sizes (min. - max.):", min(object$session$sample_size), "-", max(object$session$sample_size))
      cat("\n", "- Time points (min. - max.):", min(object$session$time_points), "-", max(object$session$time_points))

      cat("\n\n", "Suggestions:", sep = "")
      cat("\n", "- Specify the `parameter` argument in the powRICLPM_summary() function to obtain detailed results for a specific parameter.")
      cat("\n", "- Use powRICLPM_plot() to visualize the power across all conditions.")
    } else {

      # Print summary title
      cat("\n", "Post hoc powRICLPM-analysis results for parameter " , parameter, ":\n\n", sep = "")

      # Print MCMC power analysis results for specified parameter
      print(
        purrr::map_dfr(object$conditions, function(condition) {
          data.frame(sample_size = condition$sample_size,
                     time_points = condition$time_points,
                     ICC = condition$ICC,
                     condition$results[condition$results$par == parameter, -1]
          )
        })
      )
    }
  }
}

#' @title
#' Get Parameter Names From "powriclpm-Object
#'
#' @description
#' \code{\link{powRICLPM_names}} gets the names of the variables that are internally created by the powRICLPM package. Details about the naming conventions can be found in the "Details" section of \code{\link{powRICLPM}}.
#'
#' @param x A powRICLPM object.
#' @param max_set A logical indicating if the parameter names from the condition with the largest number of parameters should be returned. When simulating the power for conditions with varying number of time points, then there are different numbers of parameters over the conditions. By default, this function returns the parameter names from the condition with the smallest number of parameters, such that the returned parameter names are valid for each condition.
#'
#' @return A character vector with the names of the variables internally created by the powRICLPM package.
#'
#' @seealso \code{\link{powRICLPM_summary}}
#'
#' @examples
#' # Define population parameters for lagged effects and within-component correlations
#' Phi <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#' wSigma <- matrix(c(1, .3, .3, 1), ncol = 2, byrow = TRUE)
#'
#' # Create powRICLPM object for à priori power analysis
#' output <- powRICLPM(target_power = 0.5,
#'                      parameter = "wB2~wA1",
#'                      search_lower = 300,
#'                      search_upper = 600,
#'                      search_step = 50,
#'                      time_points = 3,
#'                      ICC = 0.5,
#'                      RI_cor = 0.3,
#'                      Phi = Phi,
#'                      wSigma = wSigma,
#'                      reps = 50,
#'                      cores = 1,
#'                      seed = 123456)
#'
#' # Get names of internally created parameters
#' powRICLPM_names(output)
#' @export
powRICLPM_names <- function(x, max_set = FALSE) {

  # Determine number of parameters per condition
  condition_length <- purrr::map_int(x$conditions, function(x) {
    length(x$results$par)
  })

  if (max_set) {

    # Return names of parameters in condition with most parameters
    return(x$conditions[[which.max(condition_length)]]$results$par)

  } else {

    # Return names of parameters in condition with least parameters
    return(x$conditions[[which.min(condition_length)]]$results$par)
  }
}
