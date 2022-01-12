#' @title
#' Create Lavaan Model Syntax for RICLPM Power Analysis
#'
#' @description
#' Creates a \pkg{lavaan} parameter table or model syntax for the specified model.
#'
#' @inheritParams powRICLPM
#' @param Psi Variance-covariance matrix of within-unit residuals from wave 2 onwards.
#' @param syntax Logical indicating whether model syntax should be created.
#'
#' @return A data frame containing the model parameters (parameter elements as characters).
#'
#' @details
#' \subsection{Data generation}
#' If lavaan model syntax needs to be created for data generation, the user must provide values for the \code{ICC}, \code{RI_cor}, \code{Phi}, \code{wSigma}, and \code{Psi} arguments. By default, these arguments are set to \code{NULL}, such that when the model syntax is made from the parameter table, these parameters are estimated rather than set.
#'
#' \subsection{Naming conventions}
#' Details on the naming conventions can be found in the "Details" section of \code{\link{powRICLPM}}.
#'
#' @examples
#' # Define population parameters for lagged effects and within-component correlations
#' Phi <- matrix(c(.5, .1, .4, .5), ncol = 2, byrow = T)
#' wSigma <- matrix(c(1 , .3, .3, 1) , ncol = 2, byrow = T)
#' Psi <- compute_Psi(Phi = Phi, wSigma = wSigma)
#'
#' # Create lavaan model syntax for data generation
#' syntax <- create_lavaan(time_points = 3,
#'                         ICC = 0.5,
#'                         RI_cor = 0.3,
#'                         Phi = Phi,
#'                         wSigma = wSigma,
#'                         Psi = Psi,
#'                         syntax = F)
create_lavaan <- function(time_points,
                          ICC = NULL,
                          RI_cor = NULL,
                          Phi = NULL,
                          wSigma = NULL,
                          Psi = NULL,
                          syntax = F) {
  # Save all arguments to list
  input <- as.list(environment())

  # Number of variables
  input$k <- 2

  # Generate default variable names
  name_var <- LETTERS[1:input$k]

  # Create matrix of names for observed variable, within, and between components
  input$name_obs <- purrr::map_dfc(name_var, paste0, 1:time_points)
  input$name_within <- purrr::map_dfc(name_var, function(x) { paste0("w", x, 1:time_points)})
  input$name_RI <- paste0("RI_", name_var)

  # Create parameter table
  lav_table <- rbind(
    lav_RI(input = input),
    lav_RI_var(input = input),
    lav_RI_cor(input = input),
    lav_within(input = input),
    lav_lagged(input = input),
    lav_within_var1(input = input),
    lav_within_cov1(input = input),
    lav_within_var2(input = input),
    lav_within_cov2(input = input),
    lav_ME(input = input)
  )

  # Remove rownames
  rownames(lav_table) <- NULL

  # Merge parameter table into model syntax
  if (syntax) {
    lav_syntax <- paste0( # Paste over parameters
      paste0( # Paste over columns
        lav_table[,1],
        lav_table[,2],
        lav_table[,3],
        lav_table[,4],
        lav_table[,5]),
      collapse = "\n")

    return(lav_syntax)
  }
  return(lav_table)
}

lav_RI <- function(input) {
  lhs <- rep(input$name_RI, each = input$time_points)
  op <- rep("=~", times = input$k*input$time_points)
  pv <- rep("1", times = input$k*input$time_points)
  con <- rep("*", times = input$k*input$time_points)
  rhs <- c(unlist(input$name_obs))
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

lav_RI_var <- function(input) {
  lhs <- rhs <- input$name_RI
  op <- rep("~~", times = input$k)
  if (is.null(input$ICC)) {
    pv <- con <- rep("", times = input$k)
  } else {
    pv <- rep((1 / (1 - input$ICC) - 1), times = input$k)
    con <- rep("*", times = input$k)
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

lav_RI_cor <- function(input) {
  # Create combinations of random intercept factors
  combnRI <- t(combn(input$name_RI, 2))

  # Create syntax (parameter table) elements
  lhs <- combnRI[, 1]
  op <- "~~"
  if (is.null(input$RI_cor)) {
    pv <- con <- ""
  } else {
    pv <- input$RI_cor
    con <- "*"
  }
  rhs <- combnRI[, 2]
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

lav_within <- function(input) {
  lhs <- c(unlist(input$name_within))
  op <- rep("=~", times = length(input$name_within))
  pv <- rep("1", times = length(input$name_within))
  con <- "*"
  rhs <- c(unlist(input$name_obs))
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

lav_lagged <- function(input) {
  lhs <- rep(c(t(input$name_within))[-(1:input$k)], each = input$k)
  op <- rep("~", times = input$k^2)
  if (is.null(input$Phi)) {
    pv <- con <- rep("", times = input$k^2)
  } else {
    pv <- c(t(input$Phi))
    con <- "*"
  }
  rhs <- c(apply(input$name_within[-input$time_points,], 1, rep, times = input$k))
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

lav_within_var1 <- function(input) {
  lhs <- rhs <- t(input$name_within[1, ])
  op <- rep("~~", times = input$k)
  if (is.null(input$wSigma)) {
    pv <- con <- rep("", times = input$k)
  } else {
    pv <- rep("1", times = input$k)
    con <- "*"
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

lav_within_cov1 <- function(input) {
  lhs <- unlist(input$name_within[1, 1])
  rhs <- unlist(input$name_within[1, 2])
  op <- rep("~~", times = 1)
  if (is.null(input$wSigma)) {
    pv <- con <- rep("", times = 1)
  } else {
    pv <- c(wSigma[lower.tri(input$wSigma)]) # Get covariances
    con <- "*"
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

lav_within_var2 <- function(input) {
  lhs <- rhs <- c(unlist(input$name_within[-1, ]))
  op <- "~~"
  if(is.null(input$Psi)){
    pv <- con <- ""
  } else {
    pv <- rep(diag(input$Psi), each = (input$time_points - 1))
    con <- "*"
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

lav_within_cov2 <- function(input) {

  # Create combinations of later within-components
  combnWComps <- t(apply(input$name_within[-1,], 1, combn, m = 2))

  # Create syntax (parameter table) elements
  lhs <- combnWComps[, 1]
  rhs <- combnWComps[, 2]
  op <- rep("~~", times = nrow(combnWComps))
  if (is.null(input$Psi)) {
    pv <- con <- rep("", nrow(combnWComps))
  } else {
    pv <- c(input$Psi[lower.tri(input$Psi)])
    con <- "*"
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

lav_ME <- function(input) {
  lhs <- rhs <- c(unlist(input$name_obs))
  op <- rep("~~", times = length(input$name_obs))
  pv <- rep("0", times = length(input$name_obs))
  con <- "*"
  return(cbind.data.frame(lhs, op, pv, con, rhs, stringsAsFactors = F))
}

