#' scenarios function (fully optimized with Rcpp)
#'
#' @param h forecast horizon
#' @param path conditional path of observables
#' @param obs position of observable(s)
#' @param shocks position of non-driving shocks (NA if all driving)
#' @param n_draws Number of posterior draws
#' @param n_sample Number of draws to sample (<= n_draws)
#' @param n_var Number of variables
#' @param n_p Number of lags
#' @param data_ Optional matrix of data (default is Z)
#' @param g Optional matrix of non-driving shocks
#' @param Sigma_g Optional covariance matrix of non-driving shocks
#'
#' @return list of mu_eps, Sigma_eps, mu_y, Sigma_y, big_b, big_M, draws_used
#' @export
#' @useDynLib APRScenario, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppProgress


scenarios <- function(h = 3, path = NULL, obs = NULL, shocks = NULL,
                     n_draws, n_sample = n_draws, n_var, n_p, data_ = Z,
                     g = NULL, Sigma_g = NULL) {

  stopifnot(length(path) == length(obs) * h)

  tmp <- big_b_and_M(h, n_draws, n_var, n_p, data_ = data_)
  big_b <- tmp[[1]]
  big_M <- tmp[[2]]

  draws_to_use <- if (n_sample < n_draws) sample(seq_len(n_draws), n_sample) else seq_len(n_draws)
  big_b <- big_b[, , draws_to_use, drop = FALSE]
  big_M <- big_M[, , draws_to_use, drop = FALSE]
  n_draws <- n_sample

  shock_idx <- if (any(is.na(shocks))) NA_integer_ else as.integer(shocks)

  # Call C++ core with optional g and Sigma_g
  out <- if (!is.null(g) && !is.null(Sigma_g)) {
    full_scenarios_core(big_b, big_M,
                        as.integer(obs),
                        as.numeric(path),
                        shock_idx,
                        h = h,
                        n_var = n_var,
                        g_ = g,
                        Sigma_g_ = Sigma_g)
  } else {
    full_scenarios_core(big_b, big_M,
                        as.integer(obs),
                        as.numeric(path),
                        shock_idx,
                        h = h,
                        n_var = n_var)
  }

  nM <- dim(big_M)[1]
  list(
    mu_eps = abind::abind(lapply(out$mu_eps, function(x) matrix(x, ncol = 1)), along = 3),
    Sigma_eps = abind::abind(lapply(out$Sigma_eps, function(x) matrix(x, nrow = nM, ncol = nM)), along = 3),
    mu_y = abind::abind(lapply(out$mu_y, function(x) matrix(x, ncol = 1)), along = 3),
    Sigma_y = abind::abind(lapply(out$Sigma_y, function(x) matrix(x, nrow = nM, ncol = nM)), along = 3),
    big_b = big_b,
    big_M = big_M,
    draws_used = draws_to_use
  )
}
