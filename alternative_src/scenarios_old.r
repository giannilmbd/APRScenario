#' scenarios function
#' @param h forecast horizon
#' @param path conditional path of observables
#' @param obs position of observable(s)
#' @param shocks position of non-driving shocks (NA if all driving)
#' @param n_draws Number of draws
#' @param n_sample Number of draws to sample (must be <= n_draws)
#' @param n_var Number of variables
#' @param n_p Number of lags
#' @param data_ (matrix optional) The data, stacking Y over X (default is Z produced by gen_mats)
#' @returns mean and variance of innovations and forecasts (four matrices) plus big_b and big_M
#' @export
scenarios_old <- function(h=3, path=NULL, obs=NULL, shocks=NULL, n_draws, n_sample=n_draws, n_var, n_p, data_=Z) {
  tmp <- big_b_and_M(h, n_draws, n_var, n_p, data_=data_)  # <<< THIS LINE
  big_b <- tmp[[1]]
  big_M <- tmp[[2]]

  ## --- sample draws ---
  draws_to_use <- if (n_sample < n_draws) sample(seq_len(n_draws), n_sample) else seq_len(n_draws)

  big_b <- big_b[,,draws_to_use, drop=FALSE]
  big_M <- big_M[,,draws_to_use, drop=FALSE]
  n_draws <- n_sample  # update internal use

  pos_cond_vars = obs
  k_0 = length(pos_cond_vars) * h
  f = array(rep(path, n_draws), dim = c(length(pos_cond_vars)*h, 1, n_draws))

  ## driving shocks
  pos_non_driving = shocks
  k_s = length(pos_non_driving) * h

  ## Selection arrays
  C = array(rep(diag(1, n_var*h), n_draws), dim = c(n_var*h, n_var*h, n_draws))
  C_h = array(C[seq(pos_cond_vars[1], n_var*h, by=n_var),,], dim = c(h, n_var*h, n_draws))
  if (length(pos_cond_vars) > 1) {
    for (cnt in 2:length(pos_cond_vars)) {
      tmp = array(C[seq(pos_cond_vars[cnt], n_var*h, by=n_var),,], dim = c(h, n_var*h, n_draws))
      C_h = abind::abind(C_h, tmp, along=1)
    }
  }

  Omega_f = parallel::mclapply(1:n_draws, FUN = function(d) {
    (C_h[,,d] %*% t(big_M[,,d])) %*% t(C_h[,,d] %*% t(big_M[,,d]))
  }, mc.cores = parallel::detectCores()) %>% simplify2array()

  if (!any(is.na(pos_non_driving))) {
    Xi_all = array(rep(diag(1, n_var*h), n_draws), dim = c(n_var*h, n_var*h, n_draws))
    Xi = array(Xi_all[seq(pos_non_driving[1], n_var*h, by=n_var),,], dim = c(h, n_var*h, n_draws))
    if (length(pos_non_driving) > 1) {
      for (cnt in 2:length(pos_non_driving)) {
        tmp = array(Xi_all[seq(pos_non_driving[cnt], n_var*h, by=n_var),,], dim = c(h, n_var*h, n_draws))
        Xi = abind::abind(Xi, tmp, along=1)
      }
    }

    C_l = parallel::mclapply(1:n_draws, FUN = function(d) {
      Xi[,,d] %*% solve(t(big_M[,,d]))
    }, mc.cores = parallel::detectCores()) %>% simplify2array()

    C_hat = abind::abind(C_h, C_l, along=1)

    tmp = parallel::mclapply(1:n_draws, FUN = function(d) {
      tmp_b <- array(big_b[,,d], dim = c(dim(big_b)[2], 1))
      C_l[,,d] %*% tmp_b
    }, mc.cores = parallel::detectCores()) %>% simplify2array() %>%
      array(dim = c(dim(C_l)[1], dim(big_b)[1], n_draws))

    f_hat = abind::abind(f, tmp, along=1)

  } else {
    C_hat = C_h
    f_hat = f
  }

  D = parallel::mclapply(1:n_draws, FUN = function(d) {
    C_hat[,,d] %*% t(big_M[,,d])
  }, mc.cores = parallel::detectCores()) %>% simplify2array()

  D_ast = parallel::mclapply(1:n_draws, FUN = function(d) {
    MASS::ginv(D[,,d])
  }, mc.cores = parallel::detectCores()) %>% simplify2array()

  if (!any(is.na(pos_non_driving))) {
    eye_s = array(rep(diag(1, k_s), n_draws), dim = c(k_s, k_s, n_draws))
    zero_0s = array(0, dim = c(k_0, k_s, n_draws))
    zero_s0 = array(0, dim = c(k_s, k_0, n_draws))

    Omega_hat = abind::abind(
      abind::abind(Omega_f, zero_0s, along=2),
      abind::abind(zero_s0, eye_s, along=2), along=1
    )
  } else {
    Omega_hat = Omega_f
  }

  eye_nh = array(rep(diag(1, n_var*h), n_draws), dim = c(n_var*h, n_var*h, n_draws))

  ## Distribution of scenario's shocks
  mu_eps <- parallel::mclapply(1:n_draws, FUN = function(d) {
    D_ast[,,d] %*% (f_hat[,,d] - C_hat[,,d] %*% array(big_b[,,d], dim = c(dim(big_b)[2], 1)))
  }, mc.cores = parallel::detectCores()) %>% simplify2array()

  Sigma_eps <- parallel::mclapply(1:n_draws, FUN = function(d) {
    D_ast[,,d] %*% Omega_hat[,,d] %*% t(D_ast[,,d]) +
      (eye_nh[,,d] - D_ast[,,d] %*% D[,,d]) %*% t(eye_nh[,,d] - D_ast[,,d] %*% D[,,d])
  }, mc.cores = parallel::detectCores()) %>% simplify2array()

  ## Distribution of scenario's outcomes
  mu_y = parallel::mclapply(1:n_draws, FUN = function(d) {
    array(big_b[,,d], dim = c(dim(big_b)[2],1)) + t(big_M[,,d]) %*% mu_eps[,,d]
  }, mc.cores = parallel::detectCores()) %>% simplify2array()

  Sigma_y = parallel::mclapply(1:n_draws, FUN = function(d) {
    t(big_M[,,d]) %*% big_M[,,d] +
      (t(big_M[,,d]) %*% D_ast[,,d]) %*% (Omega_hat[,,d] - D[,,d] %*% t(D[,,d])) %*% (t(D_ast[,,d]) %*% big_M[,,d])
  }, mc.cores = parallel::detectCores()) %>% simplify2array()

  return(list(mu_eps, Sigma_eps, mu_y, Sigma_y, big_b, big_M,draws_to_use))
}
