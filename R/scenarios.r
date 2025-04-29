#' scenarios function (optimized with Rcpp)
#' @param h forecast horizon
#' @param path conditional path of observables
#' @param obs position of observable(s)
#' @param shocks position of non-driving shocks (NA if all driving)
#' @param n_draws Number of draws
#' @param n_sample Number of draws to sample (must be <= n_draws)
#' @param n_var Number of variables
#' @param n_p Number of lags
#' @param data_ optional matrix of Y over X (default Z)
#' @return list of mu_eps, Sigma_eps, mu_y, Sigma_y, big_b, big_M, draws_used
#' @export
scenarios <- function(h = 3, path = NULL, obs = NULL, shocks = NULL,
                      n_draws, n_sample = n_draws, n_var, n_p, data_ = Z) {
  # Get big_b and big_M from prior function
  tmp <- big_b_and_M(h, n_draws, n_var, n_p, data_ = data_)
  big_b <- tmp[[1]]
  big_M <- tmp[[2]]

  draws_to_use <- if (n_sample < n_draws) sample(seq_len(n_draws), n_sample) else seq_len(n_draws)
  big_b <- big_b[,,draws_to_use, drop=FALSE]
  big_M <- big_M[,,draws_to_use, drop=FALSE]
  n_draws <- n_sample

  pos_cond_vars <- obs
  k_0 <- length(pos_cond_vars) * h
  f <- array(rep(path, n_draws), dim = c(k_0, 1, n_draws))

  # ---- Selection Matrix C_h
  row_indices <- as.vector(sapply(pos_cond_vars, function(i) seq(i, n_var*h, by=n_var)))
  C_h <- array(0, dim = c(k_0, n_var*h, n_draws))
  eye_nh <- diag(1, n_var*h)
  for (d in 1:n_draws) {
    C_h[,,d] <- eye_nh[row_indices, ]
  }

  # ---- Handle non-driving shocks
  has_non_driving <- !any(is.na(shocks))
  if (has_non_driving) {
    k_s <- length(shocks) * h
    shock_rows <- as.vector(sapply(shocks, function(i) seq(i, n_var*h, by=n_var)))
    Xi <- array(0, dim = c(k_s, n_var*h, n_draws))
    for (d in 1:n_draws) {
      Xi[,,d] <- eye_nh[shock_rows, ]
    }

    C_l <- array(0, dim = c(k_s, n_var*h, n_draws))
    for (d in 1:n_draws) {
      C_l[,,d] <- Xi[,,d] %*% solve(t(big_M[,,d]))
    }

    C_hat <- abind::abind(C_h, C_l, along = 1)

    big_b_col <- array(big_b, dim = c(dim(big_b)[2], 1, n_draws))
    tmp <- array(0, dim = c(k_s, 1, n_draws))
    for (d in 1:n_draws) {
      tmp[,,d] <- C_l[,,d] %*% big_b_col[,,d]
    }
    f_hat <- abind::abind(f, tmp, along = 1)
  } else {
    C_hat <- C_h
    f_hat <- f
  }

  # ---- D and D_ast
  D <- array(0, dim = c(dim(C_hat)[1], dim(big_M)[2], n_draws))
  for (d in 1:n_draws) {
    D[,,d] <- C_hat[,,d] %*% t(big_M[,,d])
  }

  D_ast <- array(0, dim = dim(D))
  for (d in 1:n_draws) {
    D_ast[,,d] <- MASS::ginv(D[,,d])
  }

  # ---- Omega_hat
  Omega_f <- array(0, dim = c(k_0, k_0, n_draws))
  for (d in 1:n_draws) {
    Omega_f[,,d] <- C_h[,,d] %*% t(big_M[,,d]) %*% t(C_h[,,d] %*% t(big_M[,,d]))
  }

  if (has_non_driving) {
    zero_0s <- array(0, dim = c(k_0, k_s, n_draws))
    zero_s0 <- array(0, dim = c(k_s, k_0, n_draws))
    eye_s <- array(rep(diag(1, k_s)), dim = c(k_s, k_s, n_draws))

    Omega_hat <- array(0, dim = c(k_0 + k_s, k_0 + k_s, n_draws))
    for (d in 1:n_draws) {
      Omega_hat[,,d] <- rbind(
        cbind(Omega_f[,,d], zero_0s[,,d]),
        cbind(zero_s0[,,d], eye_s[,,d])
      )
    }
  } else {
    Omega_hat <- Omega_f
  }

  # ---- Rcpp Call
  out <- scenarios_core(
    big_b_list = lapply(1:n_draws, function(d) big_b[,,d]),
    big_M_list = lapply(1:n_draws, function(d) big_M[,,d]),
    C_hat_list = lapply(1:n_draws, function(d) C_hat[,,d]),
    D_list     = lapply(1:n_draws, function(d) D[,,d]),
    D_ast_list = lapply(1:n_draws, function(d) D_ast[,,d]),
    f_hat_list = lapply(1:n_draws, function(d) f_hat[,,d]),
    Omega_hat_list = lapply(1:n_draws, function(d) Omega_hat[,,d])
  )

  return(list(
    mu_eps = out$mu_eps,
    Sigma_eps = out$Sigma_eps,
    mu_y = out$mu_y,
    Sigma_y = out$Sigma_y,
    big_b = big_b,
    big_M = big_M,
    draws_used = draws_to_use
  ))
}
