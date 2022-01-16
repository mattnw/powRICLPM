#' @title
#' Plot powRICLPM Objects
#'
#' @description
#' Plotting method for powRICLPM objects. This function visualizes results from à priori and post hoc powRICLPM analyses.
#'
#' @param x A powRICLPM object.
#' @param y Phantom argument that has no influence on the plot produced. Included here to make the generic \code{plot()} a method for objects of class powRICLPM.
#' @param parameter A character string denoting a single variable of interest. Use \code{\link{powRICLPM_names}} to get an overview of valid parameter names for the powRICLPM object.
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
#' # Visualize results
#' powRICLPM_plot(output, parameter = "wB2~wA1")
#'
#' @export
powRICLPM_plot <- function(x, y = NULL, parameter = NULL) {

  # Check arguments
  parameter <- check_parameter_summary(parameter, x)

  # Plot for à priori powRICLPM analysis
  if(x$session$type == "apriori") {

    # Combine sample sizes and simulated power across conditions
    df_plot <- data.frame(
      sample_size = purrr::map_dbl(x$conditions, function(condition) {
        condition$sample_size
        }),
      sigs = purrr::map_dbl(x$conditions, function(condition) {
        condition$results[condition$results$par == parameter, "sig"]
        })
      )

    # Create plot
    plot_output <- ggplot2::ggplot(df_plot, ggplot2::aes(x = sample_size, y = sigs)) +
      ggplot2::geom_point(shape = 19) +
      ggplot2::geom_hline(
        yintercept = x$session$target_power,
        linetype = "dashed") +
      ggplot2::labs(
        title = "Preliminary power results across different sample sizes",
        caption = paste("Results based on", x$session$reps, "replications.")) +
      ggplot2::scale_x_continuous(
        name = "Sample size",
        breaks = seq(x$session$search_lower, x$session$search_upper, x$session$search_step)) +
      ggplot2::scale_y_continuous(
        name = "Power",
        limits = c(0, 1))

  } else if (x$session$type == "posthoc") {

    # Combine sample sizes and simulated power across conditions into a data frame
    df_plot <- data.frame(
      sample_sizes = purrr::map_dbl(x$conditions, function(condition) { condition$sample_size }),
      time_points = purrr::map_dbl(x$conditions, function(condition) { condition$time_points }),
      sigs = purrr::map_dbl(x$conditions, function(condition) { condition$results[condition$results$par == parameter, "sig"] })
    )

    # Create plot
    plot_output <- ggplot2::ggplot(
      data = df_plot,
      mapping = ggplot2::aes(x = sample_sizes,
                             y = sigs,
                             color = as.factor(time_points))) +
      ggplot2::geom_point(shape = 19) +
      ggplot2::geom_line() +
      ggplot2::labs(
        title = "Simulated power across conditions",
        caption = paste("Results based on", x$session$reps, "replications."),
        color = "Number of time points") +
      ggplot2::scale_x_continuous(
        name = "Sample size",
        breaks = unique(df_plot$sample_sizes)) +
      ggplot2::scale_y_continuous(
        name = "Power",
        limits = c(0, 1)) +
      ggplot2::theme(legend.position = "bottom")
  }

  return(plot_output)
}


