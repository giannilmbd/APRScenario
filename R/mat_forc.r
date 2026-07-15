#' mat_forc function
#'##############################################################################
#'  NB: HERE WE USE Antolin-Diaz et al notation                                #
#'  B is reduced form;                                                         #
#'  A is structural;                                                           #
#'  d is intercepts                                                            #
#'  M is reduced so that E(u*u')=Sigma=(A_0*A_0')^(-1) and M_0=A_0^(-1)*Q      #
#'  Note that the code returns conflicting notation:                           #
#'  B=>A_0^(-1)*Q and                                                          #
#'  A=>B                                                                       #
#'##############################################################################
#' The recursions (Antolin-Diaz et al. 2021, Appendix A) are implemented once,
#' in the internal single-draw kernel \code{mat_forc_draw()}; \code{mat_forc()},
#' [big_b_and_M()] and [forc_h()] are thin consumers of that kernel.
#'
#' @param h (integer) forecast horison
#' @param n_draws (integer) Number of draws
#' @param n_var (integer) Number of variables
#' @param n_p (integer) Number of lags
#' @param data_ (matrix optional) The data, stacking Y over X (data and laggs)
#'        -- columns are observations (default taken from matrices$Z)
#'        NB: this is not necessarily the same as the data used to estimate the model
#'        If run counterfactuals in previoius historical period (ie not forecast) must pass the data up to previous period relative to counterfactual
#' @param matrices Optional matrices object from gen_mats() (default taken from calling environment)
#' @param max_cores number of workers used to parallelize over posterior draws
#'        (default 1 = serial; forked processes on Unix/macOS, PSOCK cluster on Windows)
#' @returns the big_b and big_M matrices of mean and IRF
#' @examples
#' \donttest{
#' library(APRScenario)
#' data(NKdata)
#'
#' # Minimal example with a toy specification
#' spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[,2:4]), p = 1)
#' est <- bsvars::estimate(spec, S = 10)  # Use small S for fast test
#' matrices<-gen_mats(posterior = est, specification = spec)
#'
#' # Example usage for matrix forecasting
#' result <- mat_forc(h = 4, n_draws = 10, n_var = 3, n_p = 1,
#'                    matrices = matrices)
#' }
#' @export
#' @import dplyr
#'
mat_forc<-function(h=1,n_draws,n_var,n_p,data_=NULL,matrices=NULL,max_cores=1){
  # Get matrices from calling environment if not provided
  if(is.null(matrices)) {
    if(exists("matrices", envir=parent.frame())) {
      matrices <- get("matrices", envir=parent.frame())
    } else {
      stop("Please provide matrices object from gen_mats() or ensure it exists in calling environment")
    }
  }

  # Get data from matrices if not provided
  if(is.null(data_)) {
    data_ <- matrices$Z
  }

  mats <- list(B_list = matrices$B_list, M = matrices$M, intercept = matrices$intercept)
  per_draw <- apply_over_draws(n_draws = n_draws, n_cores = max_cores, parallel = "auto",
                               h = h, n_var = n_var, n_p = n_p, data_ = data_, mats = mats)

  # legacy return shape: b at horizon h only (1 x n_var x n_draws),
  # M_h as a list over horizons of n_var x n_var x n_draws arrays
  b_h <- abind::abind(
    lapply(per_draw, function(x) {
      array(x$b_h[(n_var * (h - 1) + 1):(n_var * h)], dim = c(1, n_var, 1))
    }),
    along = 3
  )
  M_h <- lapply(seq_len(h), function(i)
    abind::abind(lapply(per_draw, function(x) x$M_h[[i]]), along = 3))

  return(list(b_h,M_h))
}


#' Single-draw kernel: forecast-mean terms b (all horizons, stacked) and MA/IRF
#' matrices M_h (list over horizons, M_h[[s+1]] = s steps after impact) for one
#' posterior draw d, per the recursions of Antolin-Diaz et al. (2021), Appendix A.
#' Self-contained (base R only) so it can be shipped by value to PSOCK workers.
#' @noRd
mat_forc_draw <- function(h, n_var, n_p, data_, matrices, d) {

  K_0 <- diag(n_var)
  K_h <- vector("list", h)
  K_h[[1]] <- K_0

  if (h > 1) {
    for (i in 2:h) {
      tmp2 <- matrix(0, n_var, n_var)
      for (j in 1:(i - 1)) {
        if (j <= n_p) {
          tmp1 <- matrices$B_list[[j]][, , d]
        } else {
          tmp1 <- matrix(0, n_var, n_var)
        }
        tmp2 <- tmp2 + K_h[[i - j]] %*% tmp1
      }
      K_h[[i]] <- tmp2 + K_0
    }
  }

  M_h <- vector("list", h)
  M_h[[1]] <- matrices$M[, , d]

  if (h > 1) {
    for (i in 2:h) {
      tmp2 <- matrix(0, n_var, n_var)
      # MA recursion, Antolin-Diaz et al. (2021) Appendix A: M_i = sum_j M_{i-j} B_j
      for (j in 1:min(i - 1, n_p)) {
        tmp2 <- tmp2 + M_h[[i - j]] %*% matrices$B_list[[j]][, , d]
      }
      M_h[[i]] <- tmp2
    }
  }

  N_p_list <- vector("list", n_p)

  for (l in seq_len(n_p)) {
    tmp00 <- vector("list", h)
    tmp00[[1]] <- matrices$B_list[[l]][, , d]

    if (h > 1) {
      for (i in 2:h) {
        tmp2 <- matrix(0, n_var, n_var)
        for (j in 1:min(i - 1, n_p)) {
          tmp2 <- tmp2 + tmp00[[i - j]] %*% matrices$B_list[[j]][, , d]
        }
        if ((l + i - 1) <= n_p) {
          tmp <- matrices$B_list[[l + i - 1]][, , d]
        } else {
          tmp <- matrix(0, n_var, n_var)
        }
        tmp00[[i]] <- tmp2 + tmp
      }
    }

    N_p_list[[l]] <- tmp00
  }

  b_all <- matrix(0, n_var, h)

  for (hh in seq_len(h)) {
    b_hh <- as.numeric(matrices$intercept[, d] %*% K_h[[hh]])
    for (cnt in seq_len(n_p)) {
      y_lag <- data_[(1 + n_var * (cnt - 1)):(n_var * cnt), ncol(data_)]
      b_hh <- b_hh + as.numeric(t(y_lag) %*% N_p_list[[cnt]][[hh]])
    }
    b_all[, hh] <- b_hh
  }

  list(
    b_h = as.vector(b_all),
    M_h = M_h
  )
}


#' Run mat_forc_draw() over all posterior draws with the requested backend:
#' "none" (serial), "fork" (parallel::mclapply, Unix/macOS), "psock"
#' (parallel::parLapply, all platforms incl. Windows), or "auto" (fork where
#' available, else psock). n_cores = NULL uses physical cores minus one.
#' @noRd
apply_over_draws <- function(n_draws, n_cores, parallel, h, n_var, n_p, data_, mats) {
  if (is.null(n_cores)) {
    nc <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (is.na(nc) || nc < 1L) nc <- suppressWarnings(parallel::detectCores())
    n_cores <- max(1L, min(nc - 1L, n_draws))
  }
  n_cores <- max(1L, min(as.integer(n_cores), n_draws))

  if (n_cores == 1L) parallel <- "none"
  if (parallel == "auto") {
    parallel <- if (isTRUE(capabilities("fork"))) "fork" else "psock"
  }
  if (parallel == "fork" && !isTRUE(capabilities("fork"))) {
    warning("Forking is not available on this platform (Windows); using a PSOCK cluster instead.")
    parallel <- "psock"
  }

  switch(parallel,
    none = lapply(seq_len(n_draws), mat_forc_draw,
                  h = h, n_var = n_var, n_p = n_p, data_ = data_, matrices = mats),
    fork = parallel::mclapply(seq_len(n_draws), mat_forc_draw,
                              h = h, n_var = n_var, n_p = n_p, data_ = data_, matrices = mats,
                              mc.cores = n_cores),
    psock = {
      cl <- parallel::makeCluster(n_cores, type = "PSOCK")
      on.exit(parallel::stopCluster(cl), add = TRUE)
      # self-contained copy of the kernel: workers need neither the package
      # namespace nor clusterExport
      kernel <- mat_forc_draw
      environment(kernel) <- globalenv()
      parallel::parLapply(cl, seq_len(n_draws), kernel,
                          h = h, n_var = n_var, n_p = n_p, data_ = data_, matrices = mats)
    }
  )
}
