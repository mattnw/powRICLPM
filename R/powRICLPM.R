#' @title
#' Power Analysis for the RI-CLPM
#'
#' @description
#' \code{powRICLPM()} is used to perform a Monte Carlo power analysis for the random intercept cross-lagged panel model (RI-CLPM). It can *recommend a sample size* given a desired power level and specific parameter, and compute *performance measures* (e.g., bias, mean square error, etc.) for parameters across experimental conditions. Conditions can change in terms of sample size, number of time points, proportion of between-unit variance, skewness, and kurtosis.
#'
#' @param target_power A numeric value between 0 and 1, denoting the targeted power level.
#' @param search_lower A positive integer, denoting the lower bound of the range of sample sizes to include in the power analysis.
#' @param search_upper A positive integer, denoting the upper bound of the range of sample sizes to include in the power analysis.
#' @param search_step A positive integer, denoting the increment in sample sizes.
#' @param sample_size An integer (vector) with elements at least larger than 20 (i.e., the number of free parameters in a basic bivariate RI-CLPM with 3 repeated measures), indicating the sample sizes at which to evaluate power. This argument can be specified as an alternative to the `search_` arguments to denote specific sample sizes, rather than an entire range.
#' @param time_points An integer (vector) with elements larger than 3, indicating number of time points.
#' @param ICC A numeric value denoting the proportion of variance at the between-unit level.
#' @param RI_cor A numeric value denoting the correlation between random intercepts.
#' @param Phi A matrix of standardized autoregressive and cross-lagged effects in the population. Columns represent predictors and rows represent outcomes.
#' @param wSigma A correlation matrix for the within-unit components.
#' @param skewness A numeric (vector) denoting the skewness values for the observed variables. For more information, see \code{\link[lavaan]{simulateData}}.
#' @param kurtosis A numeric (vector) denoting the kurtosis values for the observed variables. For more information, see \code{\link[lavaan]{simulateData}}.
#' @param alpha A numeric value denoting the significance criterion. It defaults to 0.05.
#' @param reps A positive integer denoting the number of Monte Carlo replications to be used during simulations.
#' @param seed An integer of length 1. If multiple cores are used, a seed of length 1 will be used to generate a full L'Ecuyer-CMRG seed for all cores. For more information, see \code{\link[furrr]{furrr_options}}.
#' @param save_path A character string naming the directory to save any (data) files to.
#'
#'
#' @details
#' \subsection{Data generation}{Data are generated using \code{\link[lavaan]{simulateData}} from the \pkg{lavaan} package. The data generating model is an RI-CLPM with a variance of 1 for the within-unit components, such that the lagged effects (specified in the \code{Phi} argument) can be interpreted as standardized effects. This implies that the process at the within-unit level must be stationary and that the \code{Phi} matrix is subject to the stationarity constraints of a VAR(1) model (i.e., eigenvalues smaller than 1). Based on the specified \code{Phi} and \code{wSigma} matrix, the residual variances and covariance between the within-unit components are computed (see \code{\link{compute_Psi}}) for wave 2 and later.}
#'
#' \subsection{Naming conventions for observed and latent variables}{The observed variables in the RI-CLPM are given default names, namely capital letters in alphabetical order, with numbers denoting the measurement occasion. For example, for a bivariate RICLPM with 3 time points, we observe \code{A1}, \code{A2}, \code{A3}, \code{B1}, \code{B2}, and \code{B3}. Their within-components are denoted by \code{wA1}, \code{wA2}, ..., \code{wB3}, respectively. The between-components have \code{RI_} prepended to the variable name, resulting in \code{RI_A} and \code{RI_B}.
#'
#' Parameters are denoted using \pkg{lavaan} model syntax. For more information, see \url{https://lavaan.ugent.be/tutorial/syntax1.html}. For example, the random intercept variances are denoted by \code{RI_A~~RI_A} and \code{RI_B~~RI_B}, the cross-lagged effects at the first wave as \code{wB2~wA1} and \code{wA2~wB1}, and the autoregressive effects as \code{wA2~wA1} and \code{wB2~wB1}. Use \code{\link{names_powRICLPM}} to extract parameter names from the `powRICLPM` object.}
#'
#' \subsection{Data analysis}{Data is analyzed using \pkg{lavaan}. To speed up the process, power analysis across different conditions can be run in parallel (implemented using \pkg{furrr}). For more information, see \url{https://jeroendmulder.github.io/powRICLPM/parallel.html}.}
#'
#' @return
#' A list containing a "conditions" and "session" element. The "condition" element is again a list, where each element is itself a list containing the input (independent variables in a power analysis context) and output (dependent variables) of the power analysis for a specific condition. This includes:
#'
#' \itemize{
#'   \item \code{sample_size}: The sample size.
#'   \item \code{time_points}: The number of time points.
#'   \item \code{ICC}: The proportion of between-unit variance.
#'   \item \code{reps}: Number of replications.
#'   \item \code{pop_synt}: \pkg{lavaan} model syntax containing population values for data generation.
#'   \item \code{pop_tab}: \pkg{lavaan} parameter table for data generation.
#'   \item \code{est_synt}: \pkg{lavaan} model syntax for estimation.
#'   \item \code{est_tab}: \pkg{lavaan} parameter table for estimation.
#'   \item \code{save_path}: The directory (data) files were saved to.
#'   \item \code{results}: Data frame containing the power analysis results.
#'   \item \code{errors}: A logical vector denoting failed Monte Carlo replications,
#'   \item \code{not_converged}: A logical vector denoting non-converged Monte Carlo replications, and
#'   \item \code{inadmissible} A logical vector denoting Monte Carlo replications that resulted in negative variances or non-positive definite matrices.
#' }
#'
#' The "session" element is a list containing information common to all conditions, including
#' \itemize{
#'   \item \code{sample_size}: The sample sizes of interest.
#'   \item \code{time_points}: The number of time points of interest.
#'   \item \code{ICC}: The proportion of between-unit variance.
#'   \item \code{reps}: The number of Monte Carlo replication.
#'   \item \code{target_power}: The desired power level.
#' }
#'
#' @author Jeroen D. Mulder \email{j.d.mulder@@uu.nl}
#' @export
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{powRICLPM_Mplus}}: Create Mplus model syntax for RI-CLPM power analysis,
#'   \item \code{\link{plot_powRICLPM}}: Visualize results of a powRICLPM analysis, and
#'   \item \code{\link{coef_powRICLPM}}: Extract performance measures for a specific parameter, across all experimental conditions.
#' }
#'
#' @examples
#' # Define population parameters for lagged effects and within-component correlations
#' Phi <- matrix(c(.4, .1, .2, .3), ncol = 2, byrow = TRUE)
#' wSigma <- matrix(c(1, .3, .3, 1), ncol = 2, byrow = TRUE)
#'
#' # Option 1 - Get a sample size recommendation
#' output1 <- powRICLPM(target_power = 0.8,
#'                      search_lower = 500,
#'                      search_upper = 600,
#'                      search_step = 50,
#'                      time_points = 3,
#'                      ICC = 0.5,
#'                      RI_cor = 0.3,
#'                      Phi = Phi,
#'                      wSigma = wSigma,
#'                      reps = 30,
#'                      seed = 123456)
#'
#' # Option 2 - Get performance measures across 6 simulation conditions
#' \dontrun{
#' output2 <- powRICLPM(sample_size = 400,
#'                      time_points = c(3, 4, 5),
#'                      ICC = c(0.3, 0.7),
#'                      RI_cor = 0.3,
#'                      Phi = Phi,
#'                      wSigma = wSigma,
#'                      reps = 100,
#'                      seed = 123456)
#'                      }
powRICLPM <- function(target_power,
                      search_lower = NULL,
                      search_upper = NULL,
                      search_step = 20,
                      sample_size = NULL,
                      time_points,
                      ICC,
                      RI_cor,
                      Phi,
                      wSigma,
                      skewness = 0,
                      kurtosis = 0,
                      alpha = 0.05,
                      reps,
                      seed = NA,
                      save_path = NA) {

  # Check arguments
  target_power <- check_target(target_power)
  time_points <- check_T(time_points)
  ICC <- check_ICC(ICC)
  RI_cor <- check_RIcor(RI_cor)
  wSigma <- check_wSigma(wSigma)
  Phi <- check_Phi(Phi)
  skewness <- check_skewness_kurtosis(skewness)
  kurtosis <- check_skewness_kurtosis(kurtosis)
  alpha <- check_alpha(alpha)
  reps <- check_reps(reps)
  seed <- check_seed(seed)

  # Check sample size
  if (is.null(sample_size)) {

    # Check arguments
    check_search(search_lower, search_upper, search_step)

    # Get candidate sample sizes
    sample_size <- seq(search_lower, search_upper, search_step)

  } else {

    # Check argument
    sample_size <- check_N(sample_size, time_points)
  }

  # Compute Psi
  Psi <- compute_Psi(Phi, wSigma)

  # Setup analysis
  object <- setup(target_power = target_power,
                  sample_size = sample_size,
                  time_points = time_points,
                  ICC = ICC,
                  RI_cor = RI_cor,
                  Phi = Phi,
                  wSigma = wSigma,
                  Psi = Psi,
                  skewness = skewness,
                  kurtosis = kurtosis,
                  alpha = alpha,
                  reps = reps,
                  seed = seed,
                  save_path = save_path)

  # Prepare progress bar
  p <- progressr::progressor(along = object$conditions)

  # Run Monte Carlo simulation for each condition
  object$conditions <- furrr::future_map(object$conditions, run_condition, p = p,
                                         .options = furrr::furrr_options(
                                           seed = seed,
                                           scheduling = 2L # Dynamic
                                         )
  )

  # Assign "powRICLPM" class to object
  class(object) <- c("powRICLPM", class(object))

  return(object)
}


