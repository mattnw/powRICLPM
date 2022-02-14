#' @title
#' Create Mplus Syntax for RICLPM Power Analysis
#'
#' @description
#' \code{powRICLPM_Mplus()} is used to create and save \code{Mplus} model syntax for the specified RI-CLPM model for performing a Monte Carlo power analysis in Mplus.
#'
#' @inheritParams powRICLPM
#' @param Psi Variance-covariance matrix of within-unit residuals from wave 2 onwards.
#'
#' @details
#' \subsection{Syntax generation}{Mplus model syntax is created in two steps: First, the \code{MODEL POPULATION} command syntax is created in which parameters are constrained to population values, and second the \code{MODEL} command syntax is created for model estimation. Ultimately, the parameter tables are pasted together to form character vectors containing the Mplus syntax to be exported.}
#'
#' \subsection{Naming conventions}{Details on the naming conventions can be found in the "Details" section of \code{\link{powRICLPM}}.}
#'
#' @return A character containing the Mplus model syntax.
#'
#' @examples
#' # Define population parameters for lagged effects and within-component correlations
#' Phi <- matrix(c(.5, .1, .4, .5), ncol = 2, byrow = TRUE)
#' wSigma <- matrix(c(1 , .3, .3, 1) , ncol = 2, byrow = TRUE)
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
  sample_size <- check_N(sample_size, time_points)
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
  input$name_obs <- suppressMessages(
    purrr::map_dfc(name_var, paste0, 1:time_points)
  )
  input$name_within <- suppressMessages(
    purrr::map_dfc(name_var, function(x) { paste0("w", x, 1:time_points)})
  )
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
    Mplus_RI_var(input, estimation = TRUE),
    Mplus_RI_cor(input, estimation = TRUE),
    Mplus_within(input),
    Mplus_lagged(input, estimation = TRUE),
    Mplus_within_var1(input, estimation = TRUE),
    Mplus_within_cov1(input, estimation = TRUE),
    Mplus_within_var2(input, estimation = TRUE),
    Mplus_within_cov2(input, estimation = TRUE),
    Mplus_ME(input),
    stringsAsFactors = FALSE)

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

Mplus_RI <- function(input) {
  lhs <- rep(input$name_RI, each = input$time_points)
  op <- " BY "
  con <- "@1"
  rhs <- c(unlist(input$name_obs))
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

Mplus_RI_var <- function(input, estimation = FALSE) {
  lhs <- input$name_RI
  op <- rhs <- rep("", times = input$k)
  if(estimation){
    con <- ""
  } else {
    con <- rep(paste0("@", 1 / (1 - input$ICC) - 1), times = input$k)
  }
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

Mplus_RI_cor <- function(input, estimation = FALSE) {

  # Create combinations of random intercept factors
  combn_RI <- t(utils::combn(input$name_RI, 2))

  # Create syntax (parameter table) elements
  lhs <- combn_RI[, 1]
  op <- rep(" WITH ", times = nrow(combn_RI))
  if(estimation){
    con <- ""
  } else {
    con <- paste0("@", input$RI_cor)
  }
  rhs <- combn_RI[, 2]
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

Mplus_within <- function(input) {
  lhs <- c(unlist(input$name_within))
  op <- rep(" BY ", times = length(input$name_within))
  con <- rep("@1", times = length(input$name_within))
  rhs <- c(unlist(input$name_obs))
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

Mplus_lagged <- function(input, estimation = FALSE) {
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

  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

Mplus_within_var1 <- function(input, estimation = FALSE) {
  lhs <- t(input$name_within[1, ])
  op <- rhs <- ""
  if(estimation){
    con <- ""
  } else {
    con <- rep("@1", times = input$k)
  }
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

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
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

Mplus_within_var2 <- function(input, estimation = FALSE){
  lhs <- c(unlist(input$name_within[-1, ]))
  op <- rhs <- ""
  if (estimation) {
    con <- ""
  } else {
    con <- rep(paste0("@", diag(input$Psi)), each = input$time_points - 1)
  }
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

Mplus_within_cov2 <- function(input, estimation = FALSE) {
  # Create within-component combinations
  combn_within <- t(apply(input$name_within[-1,], 1, utils::combn, m = 2))

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
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

Mplus_ME <- function(input) {
  lhs <- c(unlist(input$name_obs))
  op <- rhs <- ""
  con <- rep("@0", times = length(input$name_obs))
  return(cbind.data.frame(lhs, op, rhs, con, stringsAsFactors = FALSE))
}

