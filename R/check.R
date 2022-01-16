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

check_N <- function(x) {
  if(purrr::some(x, function(z) z %% 1 != 0)) {
    stop("Elements in `sample_size` should be integers.")
  }
  if(purrr::some(x, function(z) z < 0)) {
    stop("`sample_size` should be a positive integer.")
  }
  return(x)
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
  if(!is.matrix(x)) {
    stop("`Phi` should be of type `matrix`.")
  }

  x <- check_PD(x, unit = TRUE)
  return(x)
}

check_cores <- function(x) {
  # Check number of available logical cores
  available_cores <- parallelly::availableCores(logical = TRUE)

  # Check if argument matches CPU specifications
  if (!is.null(x)) {
    if (x > available_cores) {
      stop("`cores` is greater than the number of logical cores in your CPU.")
    }
    if(x < 0) {
      stop("`cores` should be a positive integer.")
    }
    if(!is.double(x)) {
      stop("`cores` should be of type `numeric`.")
    }
    if(x %% 1 != 0) {
      stop("`cores` should be an integer.")
    }
  }
  if (is.null(x)) { # Default behavior
    x <- available_cores - 1 # Leave one out by default for other processes
    warning(paste("`n_cores` not specified. Based on your machine n_cores =", x, "is chosen."))
  }
  return(x)
}

check_seed <- function(x) {
  if (!is.null(x)) {
    if(is.numeric(x)) {
      if(x %% 1 != 0) {
        stop("`seed` should be an integer.")
      }
    } else {
      stop("`seed` should be of type `numeric`.")
    }

  } else {
    warning("Be careful, no seed was specified. The exact results might not be able to be replicated.")
  }
  return(as.integer(x))
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
  if (is.null(parameter)) {

    stop("Please provide a `parameter` argument. You can use the `powRICLPM_names` function to get an overview of the parameter names in the 'powRICLPM' object.")

  } else {
    if (parameter %in% setdiff(powRICLPM_names(object), powRICLPM_names(object, max_set = TRUE))) {
      stop("The specified parameter is not a valid parameter for all conditions.")
    }
    if (!parameter %in% powRICLPM_names(object)) {
      stop("The specified parameter is not valid. Use `powRICLPM_names` to get an overview of parameter names in the 'powRICLPM' object.")
    }
    check_parameter(parameter)
  }
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
  if(!purrr::every(list(lower, upper, step), is.double)) {
    stop("`search_lower`, `search_upper`, and `search_step` should be of type `numeric`.")
  }
  if(upper < lower) {
    stop("`search_upper` should be higher than `search_lower`.")
  }
  if(step > (upper - lower)) {
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
