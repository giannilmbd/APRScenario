#' Parameter posterior updated by the scenario conditions (Waggoner-Zha / APR Algorithm 1)
#'
#' [scenarios()] computes, for every posterior draw \eqn{\theta_j \sim p(\theta \mid Y)},
#' the distribution of the future path given \eqn{\theta_j} and the conditions. Waggoner and
#' Zha (1999) and Antolin-Diaz, Petrella and Rubio-Ramirez (2021, Algorithm 1) instead
#' target the joint posterior of parameters and future path given the data *and* the
#' conditions, in which the conditions also update \eqn{\theta}. `wz_posterior()` computes
#' that posterior from the existing draws, without re-estimating the SVAR.
#'
#' Three methods are available.
#'
#' * `"hard"`: conditions \eqn{C y_f = f} on observables. The posterior is
#'   \eqn{p(\theta \mid Y, c) \propto p(\theta \mid Y)\, N(f;\, C b_\theta,\, C \Sigma_{y,\theta} C')},
#'   obtained by importance weighting of the draws. It coincides with the stationary
#'   distribution of Algorithm 1 for hard conditions on observables.
#' * `"noisy"`: the conditions are a noisy observation \eqn{f = C y_f + \eta},
#'   \eqn{\eta \sim N(0, \Omega)}. The weight is \eqn{N(f;\, C b_\theta,\, C \Sigma_{y,\theta} C' + \Omega)}
#'   and the per-draw distribution of \eqn{y_f} is the Gaussian conditional given the
#'   observation.
#' * `"apr"`: Algorithm 1 of Antolin-Diaz, Petrella and Rubio-Ramirez (2021) for whatever
#'   conditions were passed to [scenarios()], including density conditions
#'   (`Sigma_g`, or its default \eqn{\Omega_f = D D'}) and restrictions on shocks. Step 2
#'   draws \eqn{y_f} from the scenario distribution of the current draw; step 1 draws
#'   \eqn{\theta} from \eqn{p(\theta \mid Y, y_f) \propto p(\theta \mid Y)\, p(y_f \mid \theta, Y)}
#'   by resampling the existing draws. The returned weights are the Rao-Blackwellised
#'   marginal probabilities of each draw.
#'
#' The weights are relative to the draws in `scenario`, which must be equally weighted
#' draws from \eqn{p(\theta \mid Y)}. bsvarSIGNs returns such draws: with zero or narrative
#' restrictions it resamples internally and stores its own effective sample size in
#' `posterior$posterior$ess`; read it together with the `ess` reported here.
#'
#' @param scenario output of [scenarios()].
#' @param method one of `"noisy"`, `"hard"`, `"apr"`.
#' @param h forecast horizon used in [scenarios()].
#' @param obs,path conditioned variables and their path, in the format of [scenarios()].
#'   Required for `"hard"` and `"noisy"`; ignored for `"apr"`, whose conditions are those
#'   already imposed in `scenario`.
#' @param Omega \eqn{k_0 \times k_0} covariance of the observation noise for `"noisy"`,
#'   with rows ordered by variable and then horizon (the ordering of `Sigma_g` in [scenarios()]).
#' @param n_iter,burn iterations kept and discarded by the `"apr"` sampler.
#' @param resample if `TRUE`, also return `mu_y` and `Sigma_y` for draws resampled with the
#'   posterior weights, ready for [simulate_conditional_forecasts()].
#'
#' @returns A list with
#' \describe{
#'   \item{weights}{normalised posterior weight of each draw.}
#'   \item{ess}{effective sample size of the weights, \eqn{1 / \sum_j \tilde w_j^2}.}
#'   \item{kl_param}{\eqn{\sum_j \tilde w_j \log(N \tilde w_j)}: Kullback-Leibler divergence of
#'     the updated parameter posterior from \eqn{p(\theta \mid Y)} on the draws.}
#'   \item{kl_draw}{per-draw divergence of the conditional from the unconditional forecast
#'     distribution (equation 16 of Antolin-Diaz, Petrella and Rubio-Ramirez, 2021). It is
#'     `Inf` under hard conditions, whose conditional distribution is singular.}
#'   \item{kl_joint}{divergence of the joint posterior of parameters and path from the
#'     unconditional one, \eqn{\mathrm{kl\_param} + \sum_j \tilde w_j\, \mathrm{kl\_draw}_j}
#'     (chain rule).}
#'   \item{mu_y, Sigma_y}{per-draw distribution of the path under `method`.}
#'   \item{draws_resampled, mu_y_resampled, Sigma_y_resampled}{if `resample = TRUE`.}
#'   \item{ess_step, chain}{for `"apr"`: step-1 effective sample size at each iteration and
#'     the sampled draw indices.}
#' }
#'
#' @references
#' Antolin-Diaz, J., Petrella, I. and Rubio-Ramirez, J. F. (2021). Structural scenario
#' analysis with SVARs. *Journal of Monetary Economics*, 117, 798-815.
#'
#' Waggoner, D. F. and Zha, T. (1999). Conditional forecasts in dynamic multivariate models.
#' *Review of Economics and Statistics*, 81(4), 639-651.
#'
#' @examples
#' \donttest{
#' data(NKdata)
#' spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[, 2:4]), p = 1)
#' est <- bsvars::estimate(spec, S = 200, show_progress = FALSE)
#' mats <- gen_mats(posterior = est, specification = spec)
#' path <- rep(2, 4)
#' sc <- scenarios(h = 4, path = path, obs = 1, free_shocks = NA,
#'                 posterior = est, matrices = mats)
#' wz <- wz_posterior(sc, method = "hard", h = 4, obs = 1, path = path)
#' c(ess = wz$ess, kl_param = wz$kl_param)
#' }
#' @export
wz_posterior <- function(scenario, method = c("noisy", "hard", "apr"), h,
                         obs = NULL, path = NULL, Omega = NULL,
                         n_iter = 5000, burn = 500, resample = TRUE) {
  method <- match.arg(method)
  big_M <- scenario$big_M
  nh <- dim(big_M)[1]; N <- dim(big_M)[3]
  if (nh %% h != 0) stop("dim(scenario$big_M)[1] must be a multiple of h")
  n_var <- nh / h
  b <- matrix(scenario$big_b, nrow = nh, ncol = N)

  if (method %in% c("hard", "noisy")) {
    if (is.null(obs) || is.null(path)) stop("obs and path are required for method '", method, "'")
    k0 <- length(obs) * h
    f <- if (is.null(dim(path))) as.numeric(path) else as.vector(t(path))
    if (length(f) != k0) stop("path must have length(obs) * h elements")
    if (method == "noisy") {
      if (is.null(Omega)) stop("Omega is required for method 'noisy'")
      Omega <- as.matrix(Omega)
      if (!all(dim(Omega) == k0)) stop("Omega must be ", k0, " x ", k0)
    } else {
      Omega <- matrix(0, k0, k0)
    }
    C <- wz_selection(obs, n_var, h)

    lw <- numeric(N); kl_draw <- numeric(N)
    mu_y <- array(NA_real_, c(nh, 1, N)); Sigma_y <- array(NA_real_, c(nh, nh, N))
    for (d in seq_len(N)) {
      Sy <- crossprod(big_M[, , d])
      CSy <- C %*% Sy
      S <- CSy %*% t(C) + Omega
      e <- f - as.vector(C %*% b[, d])
      R <- chol(S)
      z <- backsolve(R, e, transpose = TRUE)
      lw[d] <- -sum(log(diag(R))) - 0.5 * sum(z^2)
      K <- t(backsolve(R, backsolve(R, CSy, transpose = TRUE)))
      mu <- b[, d] + as.vector(K %*% e)
      Sig <- Sy - K %*% CSy
      Sig <- (Sig + t(Sig)) / 2
      mu_y[, 1, d] <- mu; Sigma_y[, , d] <- Sig
      kl_draw[d] <- if (method == "hard") Inf else wz_kl_gauss(mu - b[, d], Sig, Sy)
    }
    w <- exp(lw - max(lw)); w <- w / sum(w)
    out <- list(method = method, weights = w)
  } else {
    if (is.null(scenario$mu_y) || is.null(scenario$Sigma_y))
      stop("scenario must be the output of scenarios()")
    mu_y <- scenario$mu_y; Sigma_y <- scenario$Sigma_y
    prep <- wz_prepare_core(big_M)
    Fc <- array(0, c(nh, nh, N))
    for (d in seq_len(N)) Fc[, , d] <- wz_psd_sqrt(Sigma_y[, , d])
    sirg <- wz_sirg_core(b, prep$W, prep$logdet, matrix(mu_y, nh, N), Fc,
                         as.integer(n_iter), as.integer(burn), sample.int(N, 1))
    w <- as.vector(sirg$pi); w <- w / sum(w)
    kl_draw <- vapply(seq_len(N), function(d)
      wz_kl_std(scenario$mu_eps[, , d], scenario$Sigma_eps[, , d]), numeric(1))
    out <- list(method = method, weights = w, ess_step = sirg$ess, chain = sirg$chain)
  }

  pos <- w > 0
  out$ess <- 1 / sum(w^2)
  out$kl_param <- sum(w[pos] * log(N * w[pos]))
  out$kl_draw <- kl_draw
  out$kl_joint <- out$kl_param + sum(w[pos] * kl_draw[pos])
  out$mu_y <- mu_y; out$Sigma_y <- Sigma_y
  if (resample) {
    idx <- sample.int(N, N, replace = TRUE, prob = w)
    out$draws_resampled <- idx
    out$mu_y_resampled <- mu_y[, , idx, drop = FALSE]
    out$Sigma_y_resampled <- Sigma_y[, , idx, drop = FALSE]
  }
  out
}

# k0 x (n_var h) selection matrix; row (j - 1) h + t picks variable obs[j] at horizon t,
# the ordering used by full_scenarios_core().
wz_selection <- function(obs, n_var, h) {
  C <- matrix(0, length(obs) * h, n_var * h)
  for (j in seq_along(obs)) for (t in seq_len(h))
    C[(j - 1) * h + t, obs[j] + (t - 1) * n_var] <- 1
  C
}

# Symmetric square-root factor F with F F' = S, for S positive semi-definite.
wz_psd_sqrt <- function(S) {
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values))
}

# KL( N(m, S1) || N(0, S0) ), S1 possibly singular (then Inf).
wz_kl_gauss <- function(m, S1, S0) {
  ev <- eigen(S1, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-10 * max(ev)) return(Inf)
  R0 <- chol(S0)
  A <- backsolve(R0, backsolve(R0, S1, transpose = TRUE))
  z <- backsolve(R0, m, transpose = TRUE)
  0.5 * (sum(diag(A)) + sum(z^2) - length(m) + 2 * sum(log(diag(R0))) - sum(log(ev)))
}

# KL( N(mu_eps, Sigma_eps) || N(0, I) ): equation 16 of APR (2021), Inf if singular.
wz_kl_std <- function(mu, S) {
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-10 * max(ev)) return(Inf)
  0.5 * (sum(ev) + sum(mu^2) - length(mu) - sum(log(ev)))
}
