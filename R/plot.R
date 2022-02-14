#' @title
#' Plot `powRICLPM` Objects
#'
#' @description
#' \code{plot_powRICLPM} visualizes results from `powRICLPM` objects using \pkg{ggplot2}.
#'
#' @param object A `powRICLPM` object.
#' @param x A character string denoting the variable to be plotted on the x-axis: Must be "sample_size", "time_points", or "ICC".
#' @param y A character string denoting the performance measure to be plotted on the y-axis. See "Details" for a list of valid performance measures.
#' @param ... Additional options that are parsed to the mapping `aes()`.
#' @param wrap A character string denoting the variable to facet the plots by.
#' @param parameter A character string denoting a single variable of interest. Use \code{\link{names_powRICLPM}} to get an overview of parameter names in the `powRICLPM` object.
#'
#' @details
#' \subsection{Performance measures}{\code{powRICLPM()} computes several performance measures that can be plotted for each condition. These include:
#' \itemize{
#'   \item \code{avg}: Average parameter estimate.
#'   \item \code{stdDev}: Standard deviation of the estimates.
#'   \item \code{SEAvg}: Mean standard error.
#'   \item \code{mse}: Mean square error.
#'   \item \code{acc}: Accuracy, the width of the confidence interval.
#'   \item \code{cover}: The proportion of times the confidence interval captures the true population value.
#'   \item \code{pwr}: The proportion of time the **p**-value is below the significance criterion.
#'   }
#' }
#'
#' \subsection{Recommendations}{It is recommended to only include study design characteristics on the x-axis (e.g., sample size and the number of repeated measures), and to facet wrap by characteristics of the data (e.g., ICC, skewness, kurtosis). This strategy emphasizes that facets represent characteristics that influence power, but cannot be influenced by the researcher. In other words, facet wraps represent different "worlds", whereas the factors on the x-axis can be tweaked by researchers to reach the desired level of power.}
#'
#' @seealso
#' \itemize{
#'   \code{\link{coef_powRICLPM}}: Extract performance measures for a specific parameter, across all experimental conditions. This function is used internally in \code{plot_powRICLPM}.
#' }
#'
#'
#' @return
#' A `ggplot2` object.
#'
#' @export
plot_powRICLPM <- function(object,
                           x,
                           y,
                           ...,
                           wrap,
                           parameter = NULL) {

  # Argument verification
  x <- match.arg(x, c("sample_size", "time_points", "ICC"))
  y <- match.arg(y, c("pwr", "avg", "stdDev", "SEAvg", "mse", "cover", "acc"))
  parameter <- check_parameter_summary(parameter, object = object)

  # Get performance table
  d <- coef_powRICLPM(object = object, parameter = parameter)

  # Create data (ggplot2-argument)
  p <- ggplot2::ggplot(d, ggplot2::aes_string(x, y, ...))

  # Create geoms
  g <- list(
    ggplot2::geom_point(shape = 19),
    ggplot2::geom_line(),
    if (y == "pwr") {ggplot2::geom_hline(yintercept = object$session$target_power, linetype = "dashed")}
  )

  # Create facets
  f <- ggplot2::facet_wrap(wrap)

  # Create scales
  s <- list(
    if (y == "pwr") {ggplot2::scale_y_continuous(name = "Power",
                                                 limits = c(0, 1))},
    if (y == "avg") {ggplot2::scale_y_continuous(name = "Average estimate")},
    if (y == "stdDev") {ggplot2::scale_y_continuous(name = "Standard error of the estimates")},
    if (y == "SEAvg") {ggplot2::scale_y_continuous(name = "Average standard error")},
    if (y == "mse") {ggplot2::scale_y_continuous(name = "Mean Square Error (MSE)")},
    if (y == "cover") {ggplot2::scale_y_continuous(name = "Coverage rate")},
    if (x == "sample_size") {ggplot2::scale_x_continuous(name = "Sample size",
                                                         breaks = d$sample_size)},
    if (x == "time_points") {ggplot2::scale_x_continuous(name = "Number of time points",
                                                         breaks = d$time_points)},
    if (x == "ICC") {ggplot2::scale_x_continuous(name = "ICC",
                                                 breaks = d$ICC)}
    )

  # Combine data, geoms, facets, and scales
  p <- p + g + f + s
  p
}





