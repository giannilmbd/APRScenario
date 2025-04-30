#' simscen function (optimized)
#'
#' @param mu_eps mean innovation
#' @param Sigma_eps variance innovation
#' @param mu_y mean forecast
#' @param Sigma_y variance forecast
#' @param big_b history forecast
#' @param big_M IRF (innovation loading)
#' @param n_sim number of simulations
#' @param h horizon
#' @param varbls variable names
#' @param idx_sampled index of random sample to use instead of full draws (from scenarios)
#'
#' @returns conditional forecast path and distribution
#'
#' @export
#' @import dplyr

SimScen <- function(mu_eps, Sigma_eps, mu_y, Sigma_y, big_b, big_M, n_sim, h, varbls, idx_sampled = 1:dim(mu_eps)[3]) {

  n_var <- length(varbls)
  n_draws_used <- length(idx_sampled)
  k <- dim(big_b)[2]  # number of states (rows of big_b)

  total_sim <- n_sim * n_draws_used

  #--- Preallocate memory
  epsilon <- array(0, dim = c(k, total_sim))
  big_b_sim <- array(0, dim = c(k, total_sim))
  big_M_sim <- array(0, dim = c(dim(big_M)[1], dim(big_M)[2], total_sim))

  #--- Draw shocks and expand big_b and big_M
  cnt <- 1
  for (i in 1:n_draws_used) {
    eps_i <- MASS::mvrnorm(n = n_sim, mu = mu_eps[,1,i], Sigma = Sigma_eps[,,i])  # (n_sim x k)
    start_idx <- (cnt-1)*n_sim + 1
    end_idx <- cnt*n_sim

    epsilon[, start_idx:end_idx] <- t(eps_i)  # (k x n_sim)
    big_b_sim[, start_idx:end_idx] <- matrix(rep(big_b[1,,i], n_sim), nrow=k)
    big_M_sim[,, start_idx:end_idx] <- array(rep(big_M[,,i], n_sim), dim = c(dim(big_M)[1], dim(big_M)[2], n_sim))

    cnt <- cnt + 1
  }

  #--- Compute y_h = big_b + M' * epsilon
  y_h <- array(0, dim = c(k, total_sim))

  for (i in 1:total_sim) {
    y_h[,i] <- big_b_sim[,i] + t(big_M_sim[,,i]) %*% epsilon[,i]
  }

  #--- Name the rows
  rownames(y_h) <- paste(rep(varbls, h), sort(rep(1:h, each=n_var)), sep = ".")

  return(y_h)
}
