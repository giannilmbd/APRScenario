# Consistency check of wz_posterior() on real bsvarSIGNs output (NKdata, the
# vignette's sign-restricted SVAR), under hard conditions on inflation.
#
# Two routes to the same posterior p(theta | Y, C y_f = f), on the same draws:
#   oracle:    method = "hard", closed-form importance weights
#              w_j ∝ N(f; C b_j, C M_j'M_j C'). Deterministic given the draws.
#   estimator: method = "apr", APR Algorithm 1 with step 1 done by resampling
#              the draws (C++ wz_sirg_core). Stochastic.
# The scenario path is the vignette's (the last h observed quarters of
# inflation), shifted by delta * s, s = sd of the unconditional forecast of
# inflation at each horizon, delta in {0, 1, 2, 3}, to see how the method
# behaves as the conditions become more informative about theta.
#
# Run from the package root, after installing the working tree into checks/lib_wz:
#   R CMD INSTALL -l checks/lib_wz . && Rscript checks/wz/wz_nkdata_check.R
# Output: checks/wz/wz_nkdata_check.log, checks/wz/wz_nkdata_check.rds
#
# TOLERANCE (fixed before the first run):
#   for each delta, over all statistics (forecast means b, reduced-form
#   coefficients A, conditional means of the unconditioned variables), the
#   z-score of (mean over R independent apr replications - hard oracle),
#   with SE = sd over replications / sqrt(R), must satisfy max |z| < 4.
#   Statistics with SE < 1e-10 must agree to 1e-8.
# Diagnostics reported, no tolerance: ESS of the oracle weights, mean step-1
# ESS of the sampler, KL_param, total-variation distance between the
# averaged sampler weights and the oracle weights, and the shift of the
# conditional forecasts relative to the at-given-theta answer.

.libPaths(c("checks/lib_wz", .libPaths()))
suppressPackageStartupMessages({
  library(APRScenario); library(bsvarSIGNs); library(bsvars); library(parallel)
})
stopifnot(packageVersion("APRScenario") >= "0.0.5.0")

ZTOL <- 4; S <- 5000; R <- 25; IT <- 2000; BURN <- 200
DELTAS <- 0:3; h <- 4; obs <- 2

log_lines <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
say("wz_nkdata_check  %s", format(Sys.time()))
say("APRScenario %s, bsvars %s, bsvarSIGNs %s", packageVersion("APRScenario"),
    packageVersion("bsvars"), packageVersion("bsvarSIGNs"))
say("S = %d posterior draws, R = %d apr replications of %d iterations (+%d burn-in)",
    S, R, IT, BURN)

## Model: the vignette's specification.
data(NKdata)
X0 <- as.data.frame(NKdata[, -1])
sr <- matrix(c(1, -1, -1,
               1,  1, -1,
               1, NA,  1), nrow = 3, byrow = TRUE)
X <- X0[1:(nrow(X0) - 4), ]
set.seed(20261120)
spec <- specify_bsvarSIGN$new(data = as.matrix(X), p = 4, sign_irf = sr)
t0 <- Sys.time()
post <- estimate(spec, S = S, show_progress = FALSE)
say("estimation: %.1f s; bsvarSIGNs ess field: %s",
    as.numeric(difftime(Sys.time(), t0, units = "secs")),
    paste(format(post$posterior$ess), collapse = " "))
mats <- gen_mats(posterior = post, specification = spec)
N <- dim(post$posterior$B)[3]; n_var <- 3

TT <- nrow(X0)
path0 <- X0[(TT - h + 1):TT, obs]
k0 <- h

## Unconditional forecast sd of inflation, from the scenario object at delta 0.
sc0 <- scenarios(h = h, path = path0, obs = obs, free_shocks = NA,
                 data_ = mats$Z, Sigma_g = matrix(0, k0, k0),
                 posterior = post, matrices = mats, n_cores = 16)
pos <- obs + n_var * (0:(h - 1))
b0 <- matrix(sc0$big_b, nrow = n_var * h)
v_unc <- rowMeans(apply(sc0$big_M, 3, function(M) diag(crossprod(M))[pos])) +
  apply(b0[pos, , drop = FALSE], 1, var)
s_unc <- sqrt(v_unc)
say("vignette path (inflation): %s", paste(sprintf("%.3f", path0), collapse = " "))
say("posterior mean forecast:   %s", paste(sprintf("%.3f", rowMeans(b0[pos, ])), collapse = " "))
say("unconditional forecast sd: %s", paste(sprintf("%.3f", s_unc), collapse = " "))

## Statistics: b (all), A (reduced form), conditional means of free elements.
A_mat <- matrix(post$posterior$A, ncol = N)
free <- setdiff(seq_len(n_var * h), pos)

RNGkind("L'Ecuyer-CMRG")
results <- list()
all_pass <- TRUE
for (delta in DELTAS) {
  f <- path0 + delta * s_unc
  sc <- if (delta == 0) sc0 else
    scenarios(h = h, path = f, obs = obs, free_shocks = NA,
              data_ = mats$Z, Sigma_g = matrix(0, k0, k0),
              posterior = post, matrices = mats, n_cores = 16)
  b <- matrix(sc$big_b, nrow = n_var * h)
  mu <- matrix(sc$mu_y, nrow = n_var * h)
  G <- rbind(b, A_mat, mu[free, , drop = FALSE])
  lab <- c(paste0("b", seq_len(nrow(b))), paste0("A", seq_len(nrow(A_mat))),
           paste0("mu_y", free))

  hard <- wz_posterior(sc, method = "hard", h = h, obs = obs, path = f, resample = FALSE)
  est_hard <- as.vector(G %*% hard$weights)

  set.seed(1000 + delta)
  t0 <- Sys.time()
  reps <- mclapply(seq_len(R), function(r)
    wz_posterior(sc, method = "apr", h = h, n_iter = IT, burn = BURN, resample = FALSE),
    mc.cores = min(R, 25), mc.set.seed = TRUE)
  t_apr <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  W <- sapply(reps, `[[`, "weights")
  E <- G %*% W
  m <- rowMeans(E); se <- apply(E, 1, sd) / sqrt(R)
  d <- m - est_hard
  z <- ifelse(se < 1e-10, ifelse(abs(d) < 1e-8, 0, Inf), d / se)
  worst <- which.max(abs(z))
  pass <- max(abs(z)) < ZTOL
  all_pass <- all_pass && pass

  # Non-vacuity: distance of the oracle from the unweighted p(theta | Y) means.
  est_unw <- rowMeans(G)
  se_unw <- apply(G, 1, sd) / sqrt(hard$ess)
  z_shift <- max(abs((est_hard - est_unw) / pmax(se_unw, 1e-12))[-(seq_len(nrow(b)))][
    seq_len(nrow(A_mat))])
  tv <- 0.5 * sum(abs(rowMeans(W) - hard$weights))
  step_ess <- mean(sapply(reps, function(x) mean(x$ess_step)))

  # Conditional forecasts of the free variables: updated vs at-given-theta.
  shift_free <- as.vector(mu[free, ] %*% hard$weights) - rowMeans(mu[free, ])

  say("")
  say("delta = %d  path = %s", delta, paste(sprintf("%.3f", f), collapse = " "))
  say("  oracle ESS %.1f of %d (%.1f%%); KL_param %.4f", hard$ess, N, 100 * hard$ess / N, hard$kl_param)
  say("  apr: mean step-1 ESS %.1f; mean ESS of averaged weights %.1f; TV(avg apr, oracle) %.4f; %.1f s",
      step_ess, 1 / sum(rowMeans(W)^2), tv, t_apr)
  say("  max |z| over %d statistics = %.2f (%s): %s", length(z), max(abs(z)), lab[worst],
      if (pass) "PASS" else "FAIL")
  say("  non-vacuity: max z of oracle vs unweighted A means = %.1f", z_shift)
  say("  shift of conditional means of free variables (updated - given theta): max |.| %.4f",
      max(abs(shift_free)))
  say("    GDP h1..h4: %s", paste(sprintf("%+.4f", shift_free[seq(1, length(free), 2)]), collapse = " "))
  say("    i   h1..h4: %s", paste(sprintf("%+.4f", shift_free[seq(2, length(free), 2)]), collapse = " "))
  results[[as.character(delta)]] <- list(f = f, ess = hard$ess, kl = hard$kl_param,
    step_ess = step_ess, tv = tv, z = z, lab = lab, pass = pass, shift_free = shift_free)
}
say("")
say("OVERALL: %s", if (all_pass) "PASS" else "FAIL")
saveRDS(results, "checks/wz/wz_nkdata_check.rds")
writeLines(log_lines, "checks/wz/wz_nkdata_check.log")
