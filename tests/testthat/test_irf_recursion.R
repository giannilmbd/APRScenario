# test_irf_recursion.R
# Guards the MA/IRF recursion (Antolin-Diaz et al. 2021, Appendix A):
# M_i = sum_{j=1}^{min(i,p)} M_{i-j} B_j.
# A mistranscription of this recursion (M_j B_j) shipped in versions <= 0.0.3.1.

test_that("mat_forc_draw reproduces brute-force VAR forecasts and IRFs", {
  set.seed(7)
  nv <- 3; np <- 2; hh <- 6; nd <- 2

  # random non-commuting lag matrices, row convention (B_list[[j]] = t(B_j))
  Bt   <- lapply(1:np, function(j) array(rnorm(nv^2 * nd, sd = 0.25), c(nv, nv, nd)))
  M0   <- array(rnorm(nv^2 * nd), c(nv, nv, nd))
  intc <- matrix(rnorm(nv * nd), nv, nd)
  yT   <- rnorm(nv); yTm1 <- rnorm(nv)
  Z    <- matrix(c(yT, yTm1), ncol = 1)      # stacked [y_T; y_{T-1}]
  mats <- list(B_list = Bt, M = M0, intercept = intc, Z = Z)

  for (d in 1:nd) {
    out <- APRScenario:::mat_forc_draw(h = hh, n_var = nv, n_p = np,
                                       data_ = Z, matrices = mats, d = d)
    b_code <- matrix(out$b_h, nv, hh)

    # brute-force mean forecast: iterate the VAR row by row
    yy <- matrix(0, np + hh, nv)
    yy[1, ] <- yTm1; yy[2, ] <- yT
    for (s in 1:hh) {
      acc <- intc[, d]
      for (j in 1:np) acc <- acc + yy[np + s - j, ] %*% Bt[[j]][, , d]
      yy[np + s, ] <- acc
    }
    expect_equal(b_code, t(yy[(np + 1):(np + hh), , drop = FALSE]),
                 tolerance = 1e-12, ignore_attr = TRUE)

    # brute-force IRFs: shock each structural innovation, iterate forward
    for (k in 1:nv) {
      irf <- matrix(0, np + hh, nv)
      eps <- rep(0, nv); eps[k] <- 1
      irf[np + 1, ] <- eps %*% M0[, , d]
      if (hh > 1) for (s in 2:hh) {
        acc <- rep(0, nv)
        for (j in 1:np) acc <- acc + irf[np + s - j, ] %*% Bt[[j]][, , d]
        irf[np + s, ] <- acc
      }
      for (s in 1:hh)
        expect_equal(out$M_h[[s]][k, ], irf[np + s, ],
                     tolerance = 1e-12, ignore_attr = TRUE)
    }
  }
})

test_that("big_b_and_M matches bsvars' impulse responses; backends agree", {
  skip_if_not_installed("bsvars")
  skip_if_not_installed("bsvarSIGNs")

  set.seed(123)
  spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[, 2:4]), p = 2)
  post <- bsvars::estimate(spec, S = 10, show_progress = FALSE)
  mats <- gen_mats(posterior = post, specification = spec)
  hh <- 5
  nd <- dim(post$posterior$B)[3]
  nv <- mats$n_var

  tmp <- big_b_and_M(h = hh, n_draws = nd, n_var = nv, n_p = mats$n_p,
                     data_ = mats$Z, matrices = mats)

  # ground truth: the authors' own IRFs on the same posterior draws;
  # big_M block (1, s+1) is the s-step response with [shock, variable] layout
  irf <- bsvars::compute_impulse_responses(post, horizon = hh - 1)
  for (s in 0:(hh - 1)) for (d in 1:nd)
    expect_equal(tmp[[2]][1:nv, (s * nv + 1):((s + 1) * nv), d],
                 t(irf[, , s + 1, d]), tolerance = 1e-8, ignore_attr = TRUE)

  # legacy per-horizon path (mat_forc) must agree with the same ground truth
  mf <- mat_forc(h = 3, n_draws = nd, n_var = nv, n_p = mats$n_p,
                 data_ = mats$Z, matrices = mats)
  for (d in 1:nd)
    expect_equal(mf[[2]][[3]][, , d], t(irf[, , 3, d]),
                 tolerance = 1e-8, ignore_attr = TRUE)

  # parallel backends reproduce the serial result exactly (2 workers: CRAN cap)
  tmp_psock <- big_b_and_M(h = hh, n_draws = nd, n_var = nv, n_p = mats$n_p,
                           data_ = mats$Z, matrices = mats,
                           n_cores = 2, parallel = "psock")
  expect_equal(tmp_psock, tmp, tolerance = 1e-12)

  if (.Platform$OS.type == "unix") {
    tmp_fork <- big_b_and_M(h = hh, n_draws = nd, n_var = nv, n_p = mats$n_p,
                            data_ = mats$Z, matrices = mats,
                            n_cores = 2, parallel = "fork")
    expect_equal(tmp_fork, tmp, tolerance = 1e-12)
  }

  # forc_h: seed-deterministic and backend-independent (2 workers: CRAN cap;
  # exercises PSOCK on Windows, fork elsewhere)
  set.seed(99)
  fc1 <- forc_h(h = 3, n_sim = 20, data_ = mats$Z, posterior = post,
                matrices = mats, max_cores = 1)
  set.seed(99)
  fc2 <- forc_h(h = 3, n_sim = 20, data_ = mats$Z, posterior = post,
                matrices = mats, max_cores = 2)
  expect_equal(fc2, fc1, tolerance = 1e-12)
})
