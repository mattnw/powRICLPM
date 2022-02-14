check_PD <- function(x, unit = FALSE){
  x_eigen <- eigen(x, only.values = TRUE)$values # Compute eigenvalues
  if (is.complex(x_eigen)){
    if (any(Re(x_eigen) < 0) | any(Im(x_eigen) < 0)) {
      stop("Complex eigenvalues: Matrix is not positive definite.")
    }
    if (unit) {
      if (any(sqrt((Re(x_eigen)^2) + (Im(x_eigen)^2) ) > 1) ) {
        stop("Complex eigenvalues are not within unit circle.")
      }
    }
  } else {
    if (any(x_eigen < 0)) {
      stop("Matrix is not positive definite.")
    }
    if (unit) {
      if (any(x_eigen >= 1)) {
        stop("Eigenvalues are not within unit circle.")
      }
    }
  }
  return(x)
}

check_K <- function(x) {
  if(!all(x %% 1 == 0)) {
    stop("Elements in `k` should be of type `numeric`.")
  }
  if (any(x < 2)) {
    stop("Elements in `k` should be larger than 1.", call. = FALSE)
  }
  return(x)
}
check_T <- function(x) {
  if(!all(x %% 1 == 0)) {
    stop("Elements in `time_points` should be integers.")
  }
  if (any(x < 3)) {
    stop("Elements in `time_points` should be larger than 2: The RI-CLPM is not identified with fewer than 3 time points.")
  }
  if (any(x > 20)) {
    warning("A large number of repeated measures can lead to computational problems. You want to consider methods for intensive longitudinal data.")
  }
  return(x)
}

check_N <- function(n, t) {
  if (purrr::some(n, function(z) z %% 1 != 0)) {
    stop("Elements in `sample_size` should be integers.")
  }
  if (purrr::some(n, function(z) z < 0)) {
    stop("`sample_size` should be a positive integer.")
  }

  k <- 2 # Number of variables (might be variable in new releases)
  t_min <- min(t) # Minimal number of repeated measures

  # Compute number of parameters
  n_parameters <- sum(factorial(1 + k) / (2*(k - 1)) * (t_min + 1),
                      k ^ 2 * (t_min - 1))

  if (purrr::some(n, function(z) z < n_parameters)) {
    stop("The number of parameters to be estimated is larger than the sample size.")
  }
  return(n)
}

check_ICC <- function(x) {
  if(!purrr::every(x, is.double)) {
    stop("Elements in `ICC` should be of type `double`.")
  }
  if(purrr::some(x, function(z) z >= 1 || z <= 0)) {
    stop("`ICC` should be between 0 and 1.")
  }
  return(x)
}

check_RIcor <- function(x) {
  if (!is.double(x)) {
    stop("`RI_cor` should be of type `double.`")
  }
  if (x < -1 || x > 1) {
    stop("`RI_cor` should be between -1 and 1 such that it can be interpreted as a correlation.")
  }
  if (length(x) > 1) {
    stop("`RI_cor` contains multiple correlations. It should be a single number.")
  }
  return(x)
}

check_wSigma <- function(x) {
  if(!is.matrix(x)) {
    stop("`wSigma` should be of type `matrix`.")
  }
  if(all(diag(x) != 1)) {
    stop("`wSigma` should be correlation matrix; 1's on the diagonal.")
  }
  x <- check_PD(x)
  return(x)
}

check_Phi <- function(x) {
  if (!is.matrix(x)) {
    stop("`Phi` should be of type `matrix`.")
  }

  x <- check_PD(x, unit = TRUE)
  return(x)
}

check_skewness_kurtosis <- function(x) {
  if (!is.numeric(x)) {
    stop("`skewness` and `kurtosis` should be of type `numeric`, `integer`, or `double`.")
  }
  return(x)
}

check_alpha <- function(x) {
  if (!is.double(x)) {
    stop("`alpha` should be of type `double`.")
  }
  if (1 < x || x < 0) {
    stop("`alpha` should be between 0 and 1.")
  }
  return(x)
}

check_seed <- function(seed) {
  if (is.na(seed)) {
    warning("No seed was specified: A seed was randomly created.")
    return(floor(stats::runif(1, -1000000, 1000000)))
  }
  if (!is.numeric(seed)) {
      stop("`seed` should be of type `numeric`.")
    } else if (seed %% 1 != 0) {
      stop("`seed` should be an integer.")
    }
  return(seed)
}

check_parameter <- function(x) {
  if (!is.null(x)) {
    if (!is.character(x)) {
      stop("`parameter` should be a character string specifying the parameter of interest.")
    }
  }
  return(x)
}

check_parameter_summary <- function(parameter, object) {
  # Include checks additional to general check_parameter()
  if (!parameter %in% names_powRICLPM(object)) {
    stop("The specified parameter is not valid. Use `names_powRICLPM()` to get an overview of parameter names in the 'powRICLPM' object.")
    }
  check_parameter(parameter)
  return(parameter)
}

check_reps <- function(x) {
  if(!is.double(x)) {
    stop("`reps` should be of type `numeric`.")
  }
  if(x %% 1 != 0 | x < 1) {
    stop("`reps` should be a positive integer.")
  }
  return(x)
}

check_search <- function(lower, upper, step) {
  if (any(purrr::map_lgl(list(lower, upper, step), is.null))) {
    stop("`search_lower`, `search_upper`, and/or `search_step` were not declared. Please provided values for all `search_` arguments.")
  }
  if (!purrr::every(list(lower, upper, step), is.double)) {
    stop("`search_lower`, `search_upper`, and `search_step` should be of type `numeric`.")
  }
  if (upper < lower) {
    stop("`search_upper` should be higher than `search_lower`.")
  }
  if (step > (upper - lower)) {
    stop("`search_step` should be smaller than or equal to the interval between `search_lower` and `search_upper`.")
  }
  return(invisible())
}

check_target <- function(x) {
  if(!is.double(x)) {
    stop("`target` should be of type `numeric`.")
  }
  if(x >= 1 | x <= 0) {
    stop("`target` should be between 0 and 1.")
  }
  return(x)
}

check_object <- function(x) {
  if (!"powRICLPM" %in% class(x)) {
    stop("`object` is not a 'powRICLPM' object.")
  }
}


