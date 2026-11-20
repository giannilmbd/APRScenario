# Validation of importance reweighting and of resampling-based APR Algorithm 1
# against an exact Waggoner-Zha / APR Gibbs sampler, in a conjugate NIW VAR(1).
#
# Run from the package root:  Rscript checks/wz/wz_conjugate_validation.R
# Output: checks/wz/wz_conjugate_validation.log, checks/wz/wz_conjugate_results.rds
#
# Methods compared for each condition type:
#   GIBBS  exact Algorithm 1: theta | Y, y_f from the NIW posterior on extended
#          data; y_f | theta from cond_moments(). The oracle.
#   IS     draws from p(theta | Y), weights p(c | theta, Y) in closed form
#          (hard and noisy only).
#   SIRG   Algorithm 1 with step 1 replaced by resampling a fixed set of draws
#          from p(theta | Y) with weights N(y_f; b_j, Sy_j) (all types).
#
# TOLERANCE (fixed before the first run): for every reported statistic,
# z = |m_method - m_gibbs| / sqrt(se_method^2 + se_gibbs^2); a method passes a
# case if max z < 4 over all statistics in that case.
# NON-VACUITY: in each case, at least one parameter's GIBBS mean must differ
# from its p(theta | Y) mean by more than 4 combined standard errors.

suppressPackageStartupMessages(library(stats))
source("checks/wz/wz_conjugate_lib.R")
set.seed(20261120)

ZTOL <- 4
N_IS <- 100000; N_SUPPORT <- 20000
N_GIBBS <- 40000; BURN_GIBBS <- 1000
N_SIRG <- 20000;  BURN_SIRG <- 500
if (nzchar(Sys.getenv("WZ_SMOKE"))) {   # smoke test only: sizes too small for the tolerance
  N_IS <- 2000; N_SUPPORT <- 1000; N_GIBBS <- 1000; BURN_GIBBS <- 100; N_SIRG <- 500; BURN_SIRG <- 50
}

## ---- data and prior --------------------------------------------------------
n <- 2; h <- 4; Tn <- 80
A_true <- matrix(c(0.7, 0.2, 0.1, 0.5), 2); c_true <- c(0.5, 0.2)
S_true <- matrix(c(1, 0.3, 0.3, 0.5), 2)
y <- matrix(0, Tn + 1, n); y[1, ] <- solve(diag(2) - A_true, c_true)
for (t in 2:(Tn + 1)) y[t, ] <- c_true + A_true %*% y[t - 1, ] + t(chol(S_true)) %*% rnorm(n)
Y <- y[-1, ]; X <- cbind(y[-(Tn + 1), ], 1); yT <- y[Tn + 1, ]
prior <- list(B0 = matrix(0, n + 1, n), V0 = diag(10, n + 1), S0 = diag(n), nu0 = n + 2)
post <- niw_post(Y, X, prior)

C <- sel_matrix(obs = 1, n = n, h = h)
theta_bar <- list(B = post$BT, Sigma = post$ST / (post$nuT - n - 1))
fm_bar <- fc_moments(theta_bar, yT, h)
sd_c <- sqrt(diag(C %*% fm_bar$Sy %*% t(C)))
f <- as.vector(C %*% fm_bar$b) + 1.5 * sd_c
Omega <- diag(0.25 * sd_c^2, length(f))

cases <- list(
  hard     = list(type = "hard",  Omega = NULL),
  noisy    = list(type = "noisy", Omega = Omega),
  apr      = list(type = "apr",   Omega = Omega),
  apr_DDt  = list(type = "apr",   Omega = NULL))

log_lines <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
say("WZ conjugate validation, %s, seed 20261120", format(Sys.time()))
say("n=%d h=%d T=%d; condition on variable 1, f = mean + 1.5 sd; Omega = 0.25 sd^2", n, h, Tn)
say("f = %s", paste(round(f, 3), collapse = " "))

## ---- check 0: forecast moments against brute-force simulation --------------
th <- niw_draw(post); fm <- fc_moments(th, yT, h); R <- 200000
sims <- matrix(0, R, n * h); yy <- matrix(yT, R, n, byrow = TRUE)
Ls <- t(chol(th$Sigma)); Am <- t(th$B[1:n, ]); cc <- th$B[n + 1, ]
for (i in 1:h) {
  yy <- sweep(yy %*% t(Am), 2, cc, "+") + matrix(rnorm(R * n), R) %*% t(Ls)
  sims[, (i - 1) * n + 1:n] <- yy
}
z_mean <- max(abs(colMeans(sims) - fm$b) / sqrt(diag(fm$Sy) / R))
cov_err <- max(abs(cov(sims) - fm$Sy)) / max(diag(fm$Sy))
say("check0 fc_moments vs simulation: max z(mean) = %.2f (tol 4); max rel cov err = %.4f (tol 0.02)",
    z_mean, cov_err)
check0_pass <- z_mean < 4 && cov_err < 0.02

## ---- unweighted posterior draws (IS proposals and SIRG support) ------------
draw_set <- function(N) {
  lapply(seq_len(N), function(i) { th <- niw_draw(post); list(th = th, fm = fc_moments(th, yT, h)) })
}
is_set  <- draw_set(N_IS)
sup_set <- draw_set(N_SUPPORT)

stats_of <- function(d, cs) theta_stats(d$th, cond_moments(d$fm, C, f, cs$type, cs$Omega))
par_names <- c(paste0("B", 1:6), "Sig11", "Sig21", "Sig22")

## SIRG precomputation: b_j and inverse upper Cholesky factors of Sy_j.
nh <- n * h
b_sup <- t(sapply(sup_set, function(d) d$fm$b))
Rinv  <- array(0, c(nh, nh, N_SUPPORT)); ldet <- numeric(N_SUPPORT)
for (j in seq_len(N_SUPPORT)) {
  Rj <- chol(sup_set[[j]]$fm$Sy); Rinv[, , j] <- backsolve(Rj, diag(nh)); ldet[j] <- sum(log(diag(Rj)))
}
loglik_y <- function(yf) {
  D <- -sweep(b_sup, 2, yf)        # rows: yf - b_j
  Z <- matrix(0, N_SUPPORT, nh)
  for (cc2 in 1:nh) { acc <- 0; for (k in 1:cc2) acc <- acc + D[, k] * Rinv[k, cc2, ]; Z[, cc2] <- acc }
  -ldet - 0.5 * rowSums(Z^2)
}

# z-score; statistics that are deterministic under both methods (e.g. the
# conditional variance of a hard-conditioned element, zero by construction)
# must agree to 1e-8 instead.
zscore <- function(m1, se1, m2, se2) {
  d <- abs(m1 - m2); s <- sqrt(se1^2 + se2^2)
  ifelse(s < 1e-10, ifelse(d < 1e-8, 0, Inf), d / pmax(s, 1e-300))
}

results <- list(); all_pass <- check0_pass
for (nm in names(cases)) {
  cs <- cases[[nm]]
  say("")
  say("=== case %s (type %s, Omega %s) ===", nm, cs$type, if (cs$type == "hard") "0" else if (is.null(cs$Omega)) "D D'" else "diag")

  # GIBBS (oracle)
  yf <- cond_moments(fm_bar, C, f, cs$type, cs$Omega)$mu
  G <- NULL; t0 <- Sys.time()
  for (it in seq_len(N_GIBBS + BURN_GIBBS)) {
    fx <- future_xy(yT, yf, n, h)
    th <- niw_draw(niw_post(rbind(Y, fx$Y), rbind(X, fx$X), prior))
    cm <- cond_moments(fc_moments(th, yT, h), C, f, cs$type, cs$Omega)
    yf <- rmvn_psd(cm$mu, cm$Sig)
    if (it > BURN_GIBBS) {
      s <- theta_stats(th, cm)
      if (is.null(G)) G <- matrix(NA_real_, N_GIBBS, length(s), dimnames = list(NULL, names(s)))
      G[it - BURN_GIBBS, ] <- s
    }
  }
  gib <- batch_mean_se(G)
  say("GIBBS done in %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))

  # Unweighted posterior (non-vacuity reference)
  G0 <- t(sapply(is_set, stats_of, cs = cs))
  unw <- list(mean = colMeans(G0), se = apply(G0, 2, sd) / sqrt(nrow(G0)))
  shift_z <- abs(gib$mean - unw$mean)[par_names] / sqrt(gib$se^2 + unw$se^2)[par_names]
  nonvac <- max(shift_z) > ZTOL
  say("non-vacuity: max z(GIBBS vs p(theta|Y)) over parameters = %.1f -> %s",
      max(shift_z), if (nonvac) "PASS" else "FAIL")

  res <- list(gibbs = gib, unweighted = unw, shift_z = shift_z, nonvacuous = nonvac)

  # IS
  if (cs$type %in% c("hard", "noisy")) {
    lw <- sapply(is_set, function(d) log_w(d$fm, C, f, cs$type, cs$Omega))
    w <- exp(lw - max(lw)); wt <- w / sum(w)
    ess <- 1 / sum(wt^2); klp <- sum(wt * log(N_IS * wt + (wt == 0)))
    isr <- is_mean_se(G0, w)
    z <- zscore(isr$mean, isr$se, gib$mean, gib$se)
    pass <- max(z) < ZTOL
    say("IS:   ESS = %.0f of %d (%.1f%%); param KL = %.3f; max z vs GIBBS = %.2f -> %s",
        ess, N_IS, 100 * ess / N_IS, klp, max(z), if (pass) "PASS" else "FAIL")
    res$is <- c(isr, list(z = z, ess = ess, kl_param = klp, pass = pass))
    all_pass <- all_pass && pass
  }

  # SIRG
  cm_sup <- lapply(sup_set, function(d) cond_moments(d$fm, C, f, cs$type, cs$Omega))
  Gs <- t(mapply(function(d, cm) theta_stats(d$th, cm), sup_set, cm_sup))
  idx <- sample.int(N_SUPPORT, 1); RB <- matrix(NA_real_, N_SIRG, ncol(Gs)); ess_it <- numeric(N_SIRG)
  t0 <- Sys.time()
  for (it in seq_len(N_SIRG + BURN_SIRG)) {
    yf <- rmvn_psd(cm_sup[[idx]]$mu, cm_sup[[idx]]$Sig)
    ll <- loglik_y(yf); w <- exp(ll - max(ll)); w <- w / sum(w)
    idx <- sample.int(N_SUPPORT, 1, prob = w)
    if (it > BURN_SIRG) { RB[it - BURN_SIRG, ] <- colSums(Gs * w); ess_it[it - BURN_SIRG] <- 1 / sum(w^2) }
  }
  colnames(RB) <- colnames(Gs)
  sirg <- batch_mean_se(RB)
  z <- zscore(sirg$mean, sirg$se, gib$mean, gib$se)
  pass <- max(z) < ZTOL
  say("SIRG: done in %.0f s; per-step ESS median %.0f (min %.0f) of %d; max z vs GIBBS = %.2f -> %s",
      as.numeric(difftime(Sys.time(), t0, units = "secs")), median(ess_it), min(ess_it),
      N_SUPPORT, max(z), if (pass) "PASS" else "FAIL")
  res$sirg <- c(sirg, list(z = z, ess_median = median(ess_it), ess_min = min(ess_it), pass = pass))
  all_pass <- all_pass && pass && nonvac

  # Effect on the scenario: forecast mean of the unconditioned variable 2.
  v2 <- paste0("mu_y", seq(2, nh, by = n))
  say("E[mu_y var2] by horizon  p(theta|Y): %s", paste(sprintf("%.3f", unw$mean[v2]), collapse = " "))
  say("                         GIBBS:      %s", paste(sprintf("%.3f", gib$mean[v2]), collapse = " "))
  results[[nm]] <- res
}

say("")
say("OVERALL: %s", if (all_pass) "PASS" else "FAIL")
writeLines(log_lines, "checks/wz/wz_conjugate_validation.log")
saveRDS(list(results = results, f = f, Omega = Omega, cases = cases), "checks/wz/wz_conjugate_results.rds")
