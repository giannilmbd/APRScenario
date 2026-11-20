# test_wz_posterior.R
# wz_posterior(): parameter posterior updated by the scenario conditions.
# Full validation against an exact Gibbs sampler: checks/wz/ (not shipped).

wz_fixture <- function(N = 60, seed = 11) {
  set.seed(seed)
  n <- 2; Tn <- 40
  Y <- matrix(rnorm(Tn * n), Tn, n)
  X <- cbind(rbind(0, Y[-Tn, ]), 1)
  A <- array(NA_real_, c(n, n + 1, N)); B <- array(NA_real_, c(n, n, N))
  for (j in seq_len(N)) {
    A[, , j] <- cbind(matrix(rnorm(n * n, sd = 0.2), n), rnorm(n, sd = 0.1))
    S <- crossprod(matrix(rnorm(n * n), n)) + diag(n)
    B[, , j] <- solve(t(chol(S)))
  }
  posterior <- list(posterior = list(A = A, B = B))
  spec <- list(data_matrices = list(Y = t(Y), X = t(X)))
  list(posterior = posterior, mats = gen_mats(posterior = posterior, specification = spec))
}

test_that("hard and noisy weights match the closed-form Gaussian likelihood", {
  fx <- wz_fixture(); h <- 3; path <- c(1, 1.5, 2)
  sc <- scenarios(h = h, path = path, obs = 1, free_shocks = NA, Sigma_g = matrix(0, h, h),
                  posterior = fx$posterior, matrices = fx$mats)
  N <- dim(sc$big_M)[3]; nh <- dim(sc$big_M)[1]
  C <- matrix(0, h, nh); for (t in 1:h) C[t, 1 + (t - 1) * 2] <- 1
  Om <- diag(0.3, h)
  for (tp in c("hard", "noisy")) {
    O <- if (tp == "noisy") Om else matrix(0, h, h)
    lw <- sapply(seq_len(N), function(d) {
      S <- C %*% crossprod(sc$big_M[, , d]) %*% t(C) + O
      e <- path - C %*% sc$big_b[1, , d]
      -0.5 * as.numeric(determinant(S)$modulus) - 0.5 * sum(e * solve(S, e))
    })
    w0 <- exp(lw - max(lw)); w0 <- w0 / sum(w0)
    wz <- wz_posterior(sc, method = tp, h = h, obs = 1, path = path,
                       Omega = if (tp == "noisy") Om else NULL)
    expect_equal(wz$weights, w0, tolerance = 1e-10)
    expect_equal(sum(wz$weights), 1)
    expect_true(wz$ess >= 1 && wz$ess <= N)
    expect_equal(wz$kl_param, sum(w0 * log(N * w0)), tolerance = 1e-10)
    expect_equal(dim(wz$mu_y_resampled), c(nh, 1, N))
  }
  # hard: per-draw conditional moments equal those of scenarios() with Sigma_g = 0
  wz <- wz_posterior(sc, method = "hard", h = h, obs = 1, path = path)
  expect_equal(wz$mu_y, sc$mu_y, tolerance = 1e-8, ignore_attr = TRUE)
  expect_true(is.infinite(wz$kl_joint))
})

test_that("apr method returns valid weights and the D D' per-draw KL", {
  fx <- wz_fixture(); h <- 3; path <- c(1, 1.5, 2)
  sc <- scenarios(h = h, path = path, obs = 1, free_shocks = NA,
                  posterior = fx$posterior, matrices = fx$mats)
  set.seed(3)
  wz <- wz_posterior(sc, method = "apr", h = h, n_iter = 300, burn = 50)
  N <- dim(sc$big_M)[3]
  expect_length(wz$weights, N)
  expect_equal(sum(wz$weights), 1)
  expect_true(all(wz$weights >= 0))
  expect_length(wz$ess_step, 300)
  expect_true(all(wz$chain >= 1 & wz$chain <= N))
  # Omega_f = D D' gives Sigma_eps = I, so eq. 16 reduces to mu_eps' mu_eps / 2
  expect_equal(wz$kl_draw, 0.5 * apply(sc$mu_eps, 3, function(x) sum(x^2)), tolerance = 1e-8)
  expect_true(is.finite(wz$kl_joint))
})
