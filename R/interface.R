#' Exported version of full_scenarios_core
#'
#' This function wraps the Rcpp-exported version of `full_scenarios_core`
#' and allows external users to call it with correct argument checks.
#'
#' @param big_b Cube of B matrices
#' @param big_M Cube of M matrices
#' @param obs Indices of constrained observables
#' @param path Flattened path for observables
#' @param shocks Indices of shocks to be recovered
#' @param h Forecast horizon
#' @param n_var Number of variables
#' @param g_ Optional vector of non-driving shocks
#' @param Sigma_g_ Optional covariance matrix of non-driving shocks
#' @export
full_scenarios_core <- function(big_b, big_M, obs, path, shocks, h, n_var, g_ = NULL, Sigma_g_ = NULL) {
  .Call(`_APRScenario_full_scenarios_core`, big_b, big_M, obs, path, shocks, h, n_var, g_, Sigma_g_)
}
