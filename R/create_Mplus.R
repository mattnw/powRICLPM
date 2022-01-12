#' @title
#' Create Mplus Syntax for RICLPM Power Analysis
#'
#' @description
#' \code{powRICLPM_Mplus()} is used to create and save \code{Mplus} model syntax for the specified RI-CLPM model for performing a Monte Carlo power analysis in Mplus.
#'
#' @inheritParams powRICLPM
#' @param Psi Variance-covariance matrix of within-unit residuals from wave 2 onwards.
#' @param syntax Logical indicating whether model syntax should be created.
#'
#' @details
#' \subsection{Syntax generation}
#' Mplus model syntax is created in two steps: First, the \code{MODEL POPULATION} command syntax is created in which parameters are constrained to population values, and second the \code{MODEL} command syntax is created for model estimation. Ultimately, the parameter tables are pasted together to form character vectors containing the Mplus syntax to be exported.
#'
#' \subsection{Naming conventions}
#' Details on the naming conventions can be found in the "Details" section of \code{\link{powRICLPM}}.
#'
#' @return A character containing the Mplus model syntax.
#'
#' @examples
#' # Define population parameters for lagged effects and within-component correlations
#' Phi <- matrix(c(.5, .1, .4, .5), ncol = 2, byrow = T)
#' wSigma <- matrix(c(1 , .3, .3, 1) , ncol = 2, byrow = T)
#' Psi <- compute_Psi(Phi = Phi, wSigma = wSigma)
#'
#' # Create and save Mplus model syntax for Monte Carlo power analysis
#' powRICLPM_Mplus(sample_size = 300,
#'                 time_points = 3,
#'                 ICC = 0.5,
#'                 RI_cor = 0.3,
#'                 Phi = Phi,
#'                 wSigma = wSigma,
#'                 reps = 10000,
#'                 seed = 123456,
#'                 save_path = "./saved")
#'
#' @export
powRICLPM_Mplus <- function(sample_size,
                            time_points,
                            ICC,
                            RI_cor,
                            Phi,
                            wSigma,
                            Psi = NULL,
                            reps = 1000,
                            seed = NULL,
                            save_path) {
  # Arguments
  sample_size <- check_N(sample_size)
  time_points <- check_T(time_points)
  ICC <- check_ICC(ICC)
  RI_cor <- check_RIcor(RI_cor)
  wSigma <- check_wSigma(wSigma)
  Phi <- check_Phi(Phi)
  reps <- check_reps(reps)
  seed <- check_seed(seed)

  # Save all arguments to list
  input <- as.list(environment())

  # Compute residual variance for within-components
  if (is.null(Psi)) {
    input$Psi <- compute_Psi(Phi = input$Phi, wSigma = input$wSigma)
  }

  # Number of variables
  input$k <- nrow(input$Phi)

  # Generate default variable names
  name_var <- LETTERS[1:input$k]

  # Create matrix of names for observed variable, within, and between components
  input$name_obs <- purrr::map_dfc(name_var, paste0, 1:time_points)
  input$name_within <- purrr::map_dfc(name_var, function(x) { paste0("w", x, 1:time_points)})
  input$name_RI <- paste0("RI_", name_var)

  # Create TITLE:, ANALYSIS:, MONTECARLO:
  preamble <- Mplus_pre(input)

  # Create MODEL POPULATION:
  Mplus_population <- rbind(
      Mplus_RI(input),
      Mplus_RI_var(input),
      Mplus_RI_cor(input),
      Mplus_within(input),
      Mplus_lagged(input),
      Mplus_within_var1(input),
      Mplus_within_cov1(input),
      Mplus_within_var2(input),
      Mplus_within_cov2(input),
      Mplus_ME(input),
      stringsAsFactors = F)

  # Create MODEL:
  Mplus_estimation <- rbind(
    Mplus_RI(input),
    Mplus_RI_var(input, estimation = T),
    Mplus_RI_cor(input, estimation = T),
    Mplus_within(input),
    Mplus_lagged(input, estimation = T),
    Mplus_within_var1(input, estimation = T),
    Mplus_within_cov1(input, estimation = T),
    Mplus_within_var2(input, estimation = T),
    Mplus_within_cov2(input, estimation = T),
    Mplus_ME(input),
    stringsAsFactors = F)

  # Add Mplus command end
  Mplus_population$end <- Mplus_estimation$end <- ";"

  # Delete Mplus_population rownames
  row.names(Mplus_population) <- row.names(Mplus_estimation) <- NULL

  # Merge parameter table into model syntax
  Mplus_syntax <- paste0(preamble, # Paste over COMMANDS
                   paste0( # Paste over parameters
                     paste0(Mplus_population[,1], # Paste over columns
                            Mplus_population[,2],
                            Mplus_population[,3],
                            Mplus_population[,4],
                            Mplus_population[,5]
                            ),
                     collapse = "\n"
                     ),
                   "\n\nMODEL:\n ",
                   paste0( # Paste over parameters
                     paste0(Mplus_estimation[,1], # Paste over columns
                            Mplus_estimation[,2],
                            Mplus_estimation[,3],
                            Mplus_estimation[,4],
                            Mplus_estimation[,5]
                            ),
                     collapse = "\n"
                     ),
                   "\n\nOUTPUT:\n TECH1 SAMPSTAT;",
                   collapse = "\n")

  # Create saving directory if it does not exist
  if (!dir.exists(save_path)) { dir.create(save_path) }

  # Save Mplus model syntax
  cat(Mplus_syntax, file = file.path(save_path, paste0("Mplus_N", sample_size, "T", time_points, ".txt")))

  # Inform user
  cat("An Mplus input file for Monte Carlo Power Analysis for the RI-CLPM has been created:")
  cat("\n - Directory: ", save_path, sep = "")
  cat("\n - File: ", paste0("Mplus_N", sample_size, "T", time_points, ".txt"), sep = "")

  # Return
  invisible()
}

#' Create Mplus Preamble for RICLPM Power Analysis
#'
#' Generate TITLE, MONTECARLO, and ANALYSIS command for RI-CLPM power analysis in Mplus.
#'
#'@param input A list with objects necessary for creating Mplus model syntax. See Details for complete list of necessary elements.
#'
#' @details
#' The \code{input} argument must containing the following elements:
#' \itemize{
#'   \item{$k}{Number of variables.}
#'   \item{$time_points}{The number of time points.}
#'   \item{$ICC}{Proportion of variance at the between-unit level.}
#'   \item{$RI_cor}{The random intercept correlations.}
#'   \item{$Phi}{Positive definite matrix of standardized autoregressive and cross-lagged effects in the population. Columns represent predictors and rows represent outcomes.}
#'   \item{wSigma}{Variance-covariance matrix of within-unit components. Variances of within-unit components must be set to 1 such that this matrix is equivalent to the correlation matrix.}
#'   \item{Psi}{Variance-covariance matrix of within-unit residuals at wave 2 and later.}
#'   \item{$name_obs}{Matrix with names of observed variables.}
#'   \item{$name_RI}{Matrix with names of random intercept factors.}
#'   \item{$name_within}{Matrix with names of within-components.}
#' }
#'
#' @return A character string containing the setup/preamble commands for power analysis in \code{Mplus}.
Mplus_pre <- function(input) {
  # Create TITLE command
  TITLE <- paste0("TITLE:\n Power analysis RICLPM with K = ", input$k, ", T = ", input$time_points, ", N = ", input$sample_size, "\n\n")

  # Create MONTECARLO command
  MONTECARLO <- paste0("MONTECARLO:\n NAMES = ", paste(unlist(input$name_obs), collapse = " "),
                       ";\n NOBSERVATIONS = ", input$sample_size,
                       ";\n NREPS = ", input$reps,
                       ";\n SEED = ", input$seed, ";\n\n")

  # Create ANALYSIS command
  ANALYSIS <- paste0("ANALYSIS:\n MODEL = NOCOV;\n\nMODEL POPULATION:\n")

  return(paste0(TITLE, MONTECARLO, ANALYSIS))
}

#' Create Mplus Syntax for Random Intercepts
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for random intercepts.
Mplus_RI <- function(input) {
  lhs <- rep(input$name_RI, each = input$time_points)
  op <- " BY "
  con <- "@1"
  rhs <- c(unlist(input$name_obs))
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

#' Create Mplus Syntax for Random Intercept Variance
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for the variances of the random intercepts.
Mplus_RI_var <- function(input, estimation = F) {
  lhs <- input$name_RI
  op <- rhs <- rep("", times = input$k)
  if(estimation){
    con <- ""
  } else {
    con <- rep(paste0("@", 1 / (1 - input$ICC) - 1), times = input$k)
  }
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

#' Create Mplus Syntax for Random Intercept Covariance
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for the covariance(s) between the random intercepts.
Mplus_RI_cor <- function(input, estimation = F) {

  # Create combinations of random intercept factors
  combn_RI <- t(combn(input$name_RI, 2))

  # Create syntax (parameter table) elements
  lhs <- combn_RI[, 1]
  op <- rep(" WITH ", times = nrow(combn_RI))
  if(estimation){
    con <- ""
  } else {
    con <- paste0("@", input$RI_cor)
  }
  rhs <- combn_RI[, 2]
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

#' Create Mplus Syntax for Within Components
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for the within components.
Mplus_within <- function(input) {
  lhs <- c(unlist(input$name_within))
  op <- rep(" BY ", times = length(input$name_within))
  con <- rep("@1", times = length(input$name_within))
  rhs <- c(unlist(input$name_obs))
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

#' Create Mplus Syntax for Within Lagged Effects
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for the within-unit lagged effects.
Mplus_lagged <- function(input, estimation = F) {
  # Create vector with outcomes
  lhs <- rep(c(t(input$name_within))[-(1:input$k)], each = input$k)

  op <- " ON "
  if(estimation){
    con <- paste0("*", rep(c(input$Phi), times = (input$time_points - 1)))
  } else {
    con <- paste0("@", rep(c(input$Phi), times = (input$time_points - 1)))
  }

  # Create vector with predictors
  rhs <- c(apply(input$name_within[-input$time_points,], 1, rep, times = input$k))

  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

#' Create Mplus Syntax for Wave 1 Within-Component Variance
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for the within-unit variance at wave 1.
Mplus_within_var1 <- function(input, estimation = F) {
  lhs <- t(input$name_within[1, ])
  op <- rhs <- ""
  if(estimation){
    con <- ""
  } else {
    con <- rep("@1", times = input$k)
  }
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

#' Create Mplus Syntax for Wave 1 Within-Component Covariance
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for the within-unit covariance at wave 1.
Mplus_within_cov1 <- function(input, estimation = F) {
  lhs <- unlist(input$name_within[1, 1])
  rhs <- unlist(input$name_within[1, 2])
  op <- rep(" WITH ", times = 1)
  if(estimation){
    con <- ""
  } else {
    resCov <- c(input$wSigma[lower.tri(input$wSigma)]) # Get covariances
    con <- paste0("@", resCov)
  }
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

#' Create Mplus Syntax for Later Residual Variance
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for within-unit residual variances at wave 2 and later.
Mplus_within_var2 <- function(input, estimation = F){
  lhs <- c(unlist(input$name_within[-1, ]))
  op <- rhs <- ""
  if (estimation) {
    con <- ""
  } else {
    con <- rep(paste0("@", diag(input$Psi)), each = input$time_points - 1)
  }
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

#' Create Mplus Syntax for Later Residual Covariance
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for the within-unit residual covariance(s) at wave 2 and later.
Mplus_within_cov2 <- function(input, estimation = F) {
  # Create within-component combinations
  combn_within <- t(apply(input$name_within[-1,], 1, combn, m = 2))

  # Create syntax (parameter table) elements
  lhs <- combn_within[, 1]
  rhs <- combn_within[, 2]
  op <- rep(" WITH ", times = nrow(combn_within))
  if (estimation) {
    con <- ""
  } else {
    resCov <- c(input$Psi[lower.tri(input$Psi)])
    con <- paste0("@", resCov)
  }
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

#' Create Mplus Syntax No Measurement Errors
#'
#' @inheritParams Mplus_pre
#'
#' @return A data frame (parameter table) with Mplus syntax for constraining the unique factor variances of the observed variables (e.g., measurement error variances) to 0.
Mplus_ME <- function(input) {
  lhs <- c(unlist(input$name_obs))
  op <- rhs <- ""
  con <- rep("@0", times = length(input$name_obs))
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = F))
}

