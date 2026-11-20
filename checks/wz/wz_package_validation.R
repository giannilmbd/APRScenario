# Validation of APRScenario::wz_posterior() and of the package path that feeds
# it (gen_mats -> big_b_and_M -> scenarios), using conjugate NIW draws passed
# to gen_mats() as a bsvarSIGNs-shaped posterior object.
#
# Run from the package root, after installing the working tree into checks/lib_wz:
#   R CMD INSTALL -l checks/lib_wz . && Rscript checks/wz/wz_package_validation.R
# Output: checks/wz/wz_package_validation.log
#
# Oracles: checks/wz/wz_conjugate_lib.R (independent R code for the forecast
# and conditional moments, checked against simulation in
# wz_conjugate_validation.R), and the exact-Gibbs means saved by that script.
#
# TOLERANCES (fixed before the first run):
#   identities (moments, weights, per-draw conditional moments, KL in y vs eps
#   coordinates): max relative difference < 1e-8;
#   posterior means vs exact Gibbs: max z < 4 over all statistics.

.libPaths(c("checks/lib_wz", .libPaths()))
suppressPackageStartupMessages(library(APRScenario))
source("checks/wz/wz_conjugate_lib.R")
oracle <- readRDS("checks/wz/wz_conjugate_results.rds")

## Data, prior, conditions: the same seed and code as wz_conjugate_validation.R.
set.seed(20261120)
n <- 2; h <- 4; Tn <- 80
A_true <- matrix(c(0.7, 0.2, 0.1, 0.5), 2); c_true <- c(0.5, 0.2)
S_true <- matrix(c(1, 0.3, 0.3, 0.5), 2)
y <- matrix(0, Tn + 1, n); y[1, ] <- solve(diag(2) - A_true, c_true)
for (t in 2:(Tn + 1)) y[t, ] <- c_true + A_true %*% y[t - 1, ] + t(chol(S_true)) %*% rnorm(n)
Y <- y[-1, ]; X <- cbind(y[-(Tn + 1), ], 1); yT <- y[Tn + 1, ]
prior <- list(B0 = matrix(0, n + 1, n), V0 = diag(10, n + 1), S0 = diag(n), nu0 = n + 2)
post <- niw_post(Y, X, prior)
f <- oracle$f; Omega <- oracle$Omega
stopifnot(isTRUE(all.equal(Y[Tn, ], yT)))

TOL_ID <- 1e-8; ZTOL <- 4
N <- 20000; N_REP <- 25; IT_REP <- 1000; BURN_REP <- 200
set.seed(7)

log_lines <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
say("wz_posterior package validation, %s; N = %d draws", format(Sys.time()), N)
relmax <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
all_pass <- TRUE
check <- function(label, val, tol) {
  ok <- is.finite(val) && val < tol
  say("%-58s %.2e (tol %.0e) %s", label, val, tol, if (ok) "PASS" else "FAIL")
  all_pass <<- all_pass && ok
}

## Conjugate draws -> bsvarSIGNs-shaped posterior.
## bsvarSIGNs: A is n x (n p + 1) with the constant last; B is the structural
## matrix with u = B^{-1} eps. gen_mats sets M = solve(t(B)), and the package's
## forecast is y_f = b + M' eps, so B = solve(t(chol(Sigma))) gives M' = chol lower.
draws <- lapply(seq_len(N), function(i) niw_draw(post))
A_arr <- array(NA_real_, c(n, n + 1, N)); B_arr <- array(NA_real_, c(n, n, N))
for (j in seq_len(N)) {
  A_arr[, , j] <- t(draws[[j]]$B)
  B_arr[, , j] <- solve(t(chol(draws[[j]]$Sigma)))
}
posterior <- list(posterior = list(A = A_arr, B = B_arr))
spec <- list(data_matrices = list(Y = t(Y), X = t(X)))
mats <- gen_mats(posterior = posterior, specification = spec)

run_scen <- function(Sigma_g)
  scenarios(h = h, path = f, obs = 1, free_shocks = NA, Sigma_g = Sigma_g,
            posterior = posterior, matrices = mats)
sc_hard <- run_scen(matrix(0, h, h))
sc_apr  <- run_scen(Omega)
sc_ddt  <- run_scen(NULL)

## 1. Package forecast moments vs the oracle library.
fms <- lapply(draws, fc_moments, yT = yT, h = h)
nh <- n * h
bpk <- matrix(sc_hard$big_b, nh, N)
check("forecast means b (package vs oracle)",
      max(sapply(seq_len(N), function(j) relmax(bpk[, j], fms[[j]]$b))), TOL_ID)
check("forecast covariance M'M (package vs oracle)",
      max(sapply(seq_len(N), function(j) relmax(crossprod(sc_hard$big_M[, , j]), fms[[j]]$Sy))), TOL_ID)

## 2. Per-draw APR conditional moments (C++ core) vs the oracle y-space formula.
C <- sel_matrix(1, n, h)
for (nm in c("apr", "ddt")) {
  sc <- if (nm == "apr") sc_apr else sc_ddt
  Om <- if (nm == "apr") Omega else NULL
  cms <- lapply(fms, cond_moments, C = C, f = f, type = "apr", Omega = Om)
  check(sprintf("scenarios() mu_y vs oracle, APR %s", nm),
        max(sapply(seq_len(N), function(j) relmax(sc$mu_y[, 1, j], cms[[j]]$mu))), TOL_ID)
  check(sprintf("scenarios() Sigma_y vs oracle, APR %s", nm),
        max(sapply(seq_len(N), function(j) relmax(sc$Sigma_y[, , j], cms[[j]]$Sig))), TOL_ID)
}

## 3. wz_posterior hard / noisy vs oracle weights and moments.
G_of <- function(cms) t(sapply(seq_len(N), function(j) theta_stats(draws[[j]], cms[[j]])))
zscore <- function(m1, se1, m2, se2) {
  d <- abs(m1 - m2); s <- sqrt(se1^2 + se2^2)
  ifelse(s < 1e-10, ifelse(d < 1e-8, 0, Inf), d / pmax(s, 1e-300))
}
for (tp in c("hard", "noisy")) {
  Om <- if (tp == "noisy") Omega else NULL
  wz <- wz_posterior(sc_hard, method = tp, h = h, obs = 1, path = f, Omega = Om)
  lw <- sapply(fms, log_w, C = C, f = f, type = tp, Omega = Om)
  w0 <- exp(lw - max(lw)); w0 <- w0 / sum(w0)
  check(sprintf("weights vs oracle, %s", tp), relmax(wz$weights, w0) * N, TOL_ID)
  cms <- lapply(fms, cond_moments, C = C, f = f, type = tp, Omega = Om)
  check(sprintf("per-draw mu_y vs oracle, %s", tp),
        max(sapply(seq_len(N), function(j) relmax(wz$mu_y[, 1, j], cms[[j]]$mu))), TOL_ID)
  check(sprintf("per-draw Sigma_y vs oracle, %s", tp),
        max(sapply(seq_len(N), function(j) relmax(wz$Sigma_y[, , j], cms[[j]]$Sig))), TOL_ID)
  if (tp == "hard")
    check("hard per-draw mu_y vs scenarios(Sigma_g = 0)",
          max(abs(wz$mu_y - sc_hard$mu_y)) / max(abs(sc_hard$mu_y)), TOL_ID)
  check(sprintf("kl_param formula, %s", tp),
        abs(wz$kl_param - sum(w0 * log(N * w0))) / max(1, wz$kl_param), TOL_ID)
  if (tp == "noisy") {
    # KL in y coordinates equals APR eq. 16 in shock coordinates eps = W (y - b).
    kl_eps <- sapply(1:200, function(j) {
      W <- solve(t(sc_hard$big_M[, , j]))
      me <- W %*% (wz$mu_y[, 1, j] - bpk[, j]); Se <- W %*% wz$Sigma_y[, , j] %*% t(W)
      0.5 * (sum(diag(Se)) + sum(me^2) - nh - as.numeric(determinant(Se)$modulus))
    })
    check("kl_draw y-space vs eps-space (APR eq. 16), noisy", relmax(wz$kl_draw[1:200], kl_eps), TOL_ID)
  }
  is <- is_mean_se(G_of(cms), wz$weights)
  gib <- oracle$results[[tp]]$gibbs
  z <- zscore(is$mean, is$se, gib$mean, gib$se)
  say("%-58s ESS %.0f (%.0f%%), kl_param %.3f, kl_joint %.3f", sprintf("wz_posterior %s:", tp),
      wz$ess, 100 * wz$ess / N, wz$kl_param, wz$kl_joint)
  check(sprintf("posterior means vs exact Gibbs, %s (max z)", tp), max(z), ZTOL)
}

## 4. wz_posterior apr (Algorithm 1 on the draws) vs exact Gibbs.
for (nm in c("apr", "apr_DDt")) {
  sc <- if (nm == "apr") sc_apr else sc_ddt
  Om <- if (nm == "apr") Omega else NULL
  cms <- lapply(fms, cond_moments, C = C, f = f, type = "apr", Omega = Om)
  G <- G_of(cms)
  reps <- t(sapply(seq_len(N_REP), function(r) {
    wz <- wz_posterior(sc, method = "apr", h = h, n_iter = IT_REP, burn = BURN_REP, resample = FALSE)
    colSums(G * wz$weights)
  }))
  m <- colMeans(reps); se <- apply(reps, 2, sd) / sqrt(N_REP)
  gib <- oracle$results[[nm]]$gibbs
  z <- zscore(m, se, gib$mean, gib$se)
  wz1 <- wz_posterior(sc, method = "apr", h = h, n_iter = 5000, burn = 500, resample = FALSE)
  say("%-58s ESS %.0f, step ESS median %.0f, kl_param %.3f, kl_joint %.3f",
      sprintf("wz_posterior apr (%s), 5000 it:", nm), wz1$ess, median(wz1$ess_step), wz1$kl_param, wz1$kl_joint)
  if (nm == "apr_DDt")
    check("kl_draw = mu_eps'mu_eps / 2 under Omega = D D'",
          relmax(wz1$kl_draw, 0.5 * apply(sc$mu_eps, 3, function(x) sum(x^2))), TOL_ID)
  check(sprintf("posterior means vs exact Gibbs, %s (max z, %d reps)", nm, N_REP), max(z), ZTOL)
}

say("OVERALL: %s", if (all_pass) "PASS" else "FAIL")
writeLines(log_lines, "checks/wz/wz_package_validation.log")
