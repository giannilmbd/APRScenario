#' simscen function
#'
#' @param mu_eps mean innovation
#' @param Sigma_eps variance innovation
#' @param mu_y mean forecast
#' @param Sigma_y variance forecast
#' @param big_b history forecast
#' @param big_M IRF (innovation loading)
#' @param n_sim number of simulations
#' @param h horizon
#' @param idx_sampled index of random sample to use instead of full draws (from scenarios.r)
#'
#' @returns conditional forecast path and distribution
#'
#' @export
#'
#' @import dplyr

SimScen_old <- function(mu_eps, Sigma_eps, mu_y, Sigma_y, big_b, big_M, n_sim, h, varbls, idx_sampled=c(1:dim(mu_eps)[3])) {

  n_var <- length(varbls)
  n_draws_used <- length(idx_sampled)  # exactly the subsampled draws

  # Draw shocks using the sampled draws
  epsilon <- parallel::mclapply(1:n_draws_used, FUN = function(i) {
    MASS::mvrnorm(n = n_sim, mu = mu_eps[,1,i], Sigma = Sigma_eps[,,i])
  }, mc.cores = parallel::detectCores() - 1) %>%
    simplify2array() %>%
    aperm(c(2,1,3)) %>%
    array(dim = c(dim(big_b)[2], n_sim * n_draws_used))

  # Expand big_b across simulations
  big_b_sim <- parallel::mclapply(1:n_draws_used, FUN = function(i) {
    array(rep(big_b[1,,i], n_sim), dim = c(dim(big_b)[2], n_sim))
  }, mc.cores = parallel::detectCores() - 1) %>%
    simplify2array() %>%
    array(dim = c(dim(big_b)[2], n_sim * n_draws_used))

  # Expand big_M across simulations
  big_M_sim <- parallel::mclapply(1:n_draws_used, FUN = function(i) {
    array(rep(big_M[,,i], n_sim), dim = c(dim(big_M)[1], dim(big_M)[2], n_sim))
  }, mc.cores = parallel::detectCores() - 1) %>%
    simplify2array() %>%
    array(dim = c(dim(big_M)[1], dim(big_M)[2], n_sim * n_draws_used))

  # Compute y_h = big_b + M' * epsilon
  y_h <- parallel::mclapply(1:(n_sim * n_draws_used), FUN = function(i) {
    big_b_sim[,i] + t(big_M_sim[,,i]) %*% epsilon[,i]
  }, mc.cores = parallel::detectCores() - 1) %>%
    simplify2array() %>%
    array(dim = c(dim(big_M)[1], n_sim * n_draws_used))

  # Name the rows
  rownames(y_h) <- paste(rep(varbls, h), sort(rep(1:h, n_var)), sep = '.')

  return(y_h)
}
