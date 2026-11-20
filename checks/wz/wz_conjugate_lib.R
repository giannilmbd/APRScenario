# Conjugate Normal-inverse-Wishart VAR(p=1) with constant: building blocks for
# the Waggoner-Zha / APR (2021) Algorithm 1 validation.
#
# Notation. Data rows y_t' = x_t' B + u_t', x_t = (y_{t-1}', 1)', u_t ~ N(0, Sigma).
# B is k x n (k = n + 1); A = t(B[1:n, ]) is the lag matrix, cc = B[n + 1, ].
# Prior: Sigma ~ IW(S0, nu0), vec(B) | Sigma ~ N(vec(B0), Sigma %x% V0).
# Future path y_f = (y_{T+1}', ..., y_{T+h}')' (nh x 1), stacked by horizon.
# Conditions on observables: C y_f with C the k0 x nh selection matrix.

niw_post <- function(Y, X, prior) {
  V0i <- solve(prior$V0)
  VT  <- solve(V0i + crossprod(X))
  BT  <- VT %*% (V0i %*% prior$B0 + crossprod(X, Y))
  ST  <- prior$S0 + crossprod(Y) + t(prior$B0) %*% V0i %*% prior$B0 -
         t(BT) %*% solve(VT) %*% BT
  list(BT = BT, VT = VT, ST = (ST + t(ST)) / 2, nuT = prior$nu0 + nrow(Y))
}

niw_draw <- function(post) {
  n <- ncol(post$BT)
  W <- stats::rWishart(1, post$nuT, solve(post$ST))[, , 1]
  Sigma <- solve(W); Sigma <- (Sigma + t(Sigma)) / 2
  Z <- matrix(stats::rnorm(length(post$BT)), nrow(post$BT), n)
  B <- post$BT + t(chol(post$VT)) %*% Z %*% chol(Sigma)
  list(B = B, Sigma = Sigma)
}

# Unconditional forecast moments of y_f given theta and y_T:
# y_f = b + Phi u_f, Sigma_y = Phi (I_h %x% Sigma) Phi'.
fc_moments <- function(theta, yT, h) {
  n <- length(yT); A <- t(theta$B[1:n, , drop = FALSE]); cc <- theta$B[n + 1, ]
  b <- numeric(n * h); y <- yT
  Apow <- vector("list", h); Apow[[1]] <- diag(n)
  for (i in seq_len(h)) {
    y <- cc + A %*% y; b[(i - 1) * n + 1:n] <- y
    if (i > 1) Apow[[i]] <- A %*% Apow[[i - 1]]
  }
  Phi <- matrix(0, n * h, n * h)
  for (i in seq_len(h)) for (j in seq_len(i))
    Phi[(i - 1) * n + 1:n, (j - 1) * n + 1:n] <- Apow[[i - j + 1]]
  Sy <- Phi %*% kronecker(diag(h), theta$Sigma) %*% t(Phi)
  list(b = b, Sy = (Sy + t(Sy)) / 2)
}

sel_matrix <- function(obs, n, h) {
  C <- matrix(0, length(obs) * h, n * h)
  r <- 0
  for (v in obs) for (t in seq_len(h)) { r <- r + 1; C[r, (t - 1) * n + v] <- 1 }
  C
}

# Conditional distribution of y_f given theta under the three condition types.
#  "hard":  C y_f = f                                    (Omega = 0)
#  "noisy": f = C y_f + eta, eta ~ N(0, Omega)           (Bayesian signal)
#  "apr":   C y_f ~ N(f, Omega), APR (2021) eqs. (13)-(14) written in y-space:
#           mu = b + K (f - C b),  Sig = Sy + K (Omega - C Sy C') K',
#           K  = Sy C' (C Sy C')^{-1}.  Omega = NULL means Omega = C Sy C' (= D D').
cond_moments <- function(fm, C, f, type, Omega = NULL) {
  CSC <- C %*% fm$Sy %*% t(C)
  if (type == "hard") {
    K <- fm$Sy %*% t(C) %*% solve(CSC)
    Sig <- fm$Sy - K %*% C %*% fm$Sy
  } else if (type == "noisy") {
    K <- fm$Sy %*% t(C) %*% solve(CSC + Omega)
    Sig <- fm$Sy - K %*% C %*% fm$Sy
  } else if (type == "apr") {
    K <- fm$Sy %*% t(C) %*% solve(CSC)
    Om <- if (is.null(Omega)) CSC else Omega
    Sig <- fm$Sy + K %*% (Om - CSC) %*% t(K)
  } else stop("unknown type")
  list(mu = as.vector(fm$b + K %*% (f - C %*% fm$b)), Sig = (Sig + t(Sig)) / 2)
}

# Draw from N(mu, Sig) with Sig possibly singular (hard conditions).
rmvn_psd <- function(mu, Sig) {
  e <- eigen(Sig, symmetric = TRUE)
  ev <- pmax(e$values, 0)
  as.vector(mu + e$vectors %*% (sqrt(ev) * stats::rnorm(length(mu))))
}

# log N(x; m, S) for nonsingular S.
ldmvn <- function(x, m, S) {
  R <- chol(S); z <- backsolve(R, x - m, transpose = TRUE)
  -0.5 * length(x) * log(2 * pi) - sum(log(diag(R))) - 0.5 * sum(z^2)
}

# Log importance weight log p(c | theta, Y).
log_w <- function(fm, C, f, type, Omega = NULL) {
  CSC <- C %*% fm$Sy %*% t(C)
  S <- if (type == "hard") CSC else if (type == "noisy") CSC + Omega else
    stop("closed-form weight exists only for hard and noisy conditions")
  ldmvn(as.vector(f), as.vector(C %*% fm$b), S)
}

# Lag-1 regressor rows for the future path (Gibbs step 1 on extended data).
future_xy <- function(yT, yf, n, h) {
  Yf <- matrix(yf, h, n, byrow = TRUE)
  lags <- rbind(yT, Yf[-h, , drop = FALSE])
  list(Y = Yf, X = cbind(lags, 1))
}

# Scalar functions of theta to compare across methods: parameters, and the
# per-draw conditional forecast mean and variance of every element of y_f.
theta_stats <- function(theta, cm) {
  n <- ncol(theta$B)
  c(setNames(as.vector(theta$B), paste0("B", seq_along(theta$B))),
    Sig11 = theta$Sigma[1, 1], Sig21 = theta$Sigma[2, 1], Sig22 = theta$Sigma[2, 2],
    setNames(cm$mu, paste0("mu_y", seq_along(cm$mu))),
    setNames(diag(cm$Sig), paste0("var_y", seq_along(cm$mu))))
}

# Monte Carlo standard errors.
is_mean_se <- function(G, w) {
  wt <- w / sum(w); m <- colSums(G * wt)
  se <- sqrt(colSums(wt^2 * sweep(G, 2, m)^2))
  list(mean = m, se = se)
}
batch_mean_se <- function(G, nb = 50) {
  N <- nrow(G); L <- floor(N / nb); G <- G[seq_len(nb * L), , drop = FALSE]
  bm <- apply(G, 2, function(g) colMeans(matrix(g, L, nb)))
  list(mean = colMeans(G), se = apply(bm, 2, sd) / sqrt(nb))
}
