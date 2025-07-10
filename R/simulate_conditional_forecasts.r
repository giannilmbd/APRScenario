#' Simulate paths from conditional forecast distributions
#'
#' @param mu_y  Array (n_state × 1 × n_draws): conditional forecast mean
#' @param Sigma_y Array (n_state × n_state × n_draws): conditional forecast variance
#' @param varnames Character vector of variable names (length = number of variables)
#' @param n_sim Number of simulations per draw
#'
#' @return Array of dimensions (n_state x n_sim x n_draws) of simulated draws with named rows
#' @note Users should set their own seed before calling this function if reproducible results are desired.
#' @examples
#' # Example with simulated data
#' # Create example data dimensions
#' n_var <- 3
#' h <- 2
#' n_draws <- 5
#' n_state <- n_var * h
#' 
#' # Simulate conditional forecast means and covariances
#' set.seed(123)
#' mu_y <- array(rnorm(n_state * 1 * n_draws), dim = c(n_state, 1, n_draws))
#' Sigma_y <- array(0, dim = c(n_state, n_state, n_draws))
#' for (d in 1:n_draws) {
#'   temp_cov <- matrix(rnorm(n_state^2), n_state, n_state)
#'   Sigma_y[,,d] <- temp_cov %*% t(temp_cov) + diag(n_state) * 0.1
#' }
#' 
#' # Variable names
#' varnames <- c("GDP", "CPI", "FFR")
#' 
#' # Simulate conditional forecasts
#' sims <- simulate_conditional_forecasts(mu_y, Sigma_y, varnames, n_sim = 50)
#' print(dim(sims))
#' print(rownames(sims)[1:6])
#' @export
simulate_conditional_forecasts <- function(mu_y, Sigma_y, varnames, n_sim = 1000) {
  stopifnot(length(dim(mu_y)) == 3, length(dim(Sigma_y)) == 3)

  n_state <- dim(mu_y)[1]
  n_draws <- dim(mu_y)[3]
  n_var <- length(varnames)
  h <- n_state / n_var
  stopifnot(h == floor(h))

  # Construct state names: var.1, var.2, ..., preserving SVAR order
  state_names <- paste0(rep(varnames, times = h), ".", rep(1:h, each = n_var))

  sims <- array(NA_real_, dim = c(n_state, n_sim, n_draws),
                dimnames = list(state_names, NULL, NULL))
  
  for (d in seq_len(n_draws)) {
    sims[,,d] <- t(MASS::mvrnorm(n_sim, mu = mu_y[,1,d], Sigma = Sigma_y[,,d]))
  }
  sim_red<-matrix(aperm(sims,c(1,3,2)),nrow=dim(sims)[1])
  rownames(sim_red)<-state_names
  return(sim_red)
}
