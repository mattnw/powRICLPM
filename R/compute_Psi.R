#' @title
#' Compute Residual Variances of Lagged Within-Components
#'
#' @description
#' Within an RI-CLPM context, this function computes the variance-covariance matrix of the within-unit residuals from wave 2 and later, given the lagged effects in \code{Phi} and an observed variance-covariance matrix \code{Sigma}.
#'
#' @inheritParams powRICLPM
#'
#' @return A variance-covariance matrix for within-unit residuals from wave 2 and later.
#'
#' @details
#' The function is based on Equation (3.26) in Kim and Nelson (1999, p. 27).
#'
#' @examples
#' # Set lagged effects
#' Phi = matrix(c(.2, .15, .10, .3), ncol = 2, byrow = T)
#'
#' # Set variance-covariance of process
#' wSigma = matrix(c(1 , .3, .3, 1) , ncol = 2, byrow = T)
#'
#' # Compute residual (co)variances
#' compute_resVar(Phi = Phi, wSigma = wSigma)
compute_Psi <- function(Phi, wSigma){
  matrix((diag(length(wSigma)) - Phi %x% Phi) %*% c(wSigma), nrow = nrow(Phi))
}
