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
#' @return A list with elements depending on the input configuration. Typically includes:
#' \describe{
#'   \item{mu_eps}{Matrix of mean structural shocks}
#'   \item{Sigma_eps}{Covariance matrix of structural shocks}
#'   \item{mu_y}{Matrix of conditional means of observables}
#'   \item{Sigma_y}{Covariance matrix of observables}
#'   \item{big_b}{Slice of B matrices used}
#'   \item{big_M}{Slice of M matrices used}
#'   \item{draws_used}{Indices of posterior draws used in the simulation}
#' }
#' @examples
#' \dontrun{
#' # This function is typically called internally by scenarios()
#' # Example usage with simulated data:
#' big_b <- array(rnorm(9*4*10), dim = c(9, 4, 10))
#' big_M <- array(rnorm(9*9*10), dim = c(9, 9, 10))
#' result <- full_scenarios_core(big_b, big_M, obs = 1:2, 
#'                               path = c(1.0, 1.1), shocks = NA, 
#'                               h = 2, n_var = 3)
#' }
#' @export
full_scenarios_core <- function(big_b, big_M, obs, path, shocks, h, n_var, g_ = NULL, Sigma_g_ = NULL) {
  .Call(`_APRScenario_full_scenarios_core`, big_b, big_M, obs, path, shocks, h, n_var, g_, Sigma_g_)
}
