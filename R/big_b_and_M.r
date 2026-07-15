#' big_b_and_M
#' This function returns the extended b and M matrices as in APR
#'
#' Computation is done draw by draw (all horizons in a single pass) and can be
#' parallelized over posterior draws on every platform: forked processes
#' (\code{parallel::mclapply}) on Unix/macOS, a PSOCK cluster
#' (\code{parallel::parLapply}) on Windows.
#'
#' @param h forecast horison
#' @param n_draws Number of draws
#' @param n_var Number of variables
#' @param n_p Number of lags
#' @param data_ (matrix optional) The data, stacking Y over X (data and laggs)
#'        -- columns are observations (default taken from matrices$Z)
#'        NB: this is not necessarily the same as the data used to estimate the model
#'        If run counterfactuals in previoius historical period (ie not forecast) must pass the data up to previous period relative to counterfactual
#' @param matrices Optional matrices object from gen_mats() (default taken from calling environment)
#' @param n_cores Number of parallel workers (default 1 = serial). If NULL, uses
#'        the number of physical cores minus one, capped by \code{n_draws}.
#' @param parallel Parallel backend: \code{"auto"} (default; fork where available,
#'        i.e. Unix/macOS, otherwise PSOCK), \code{"fork"}, \code{"psock"}, or
#'        \code{"none"}. Ignored when \code{n_cores = 1}.
#' @returns the big_b and big_M matrices of mean and IRF
#' @examples
#' \dontrun{
#' # Example usage for creating extended matrices (serial)
#' result <- big_b_and_M(h = 4, n_draws = 1000, n_var = 3, n_p = 2,
#'                       matrices = matrices)
#' # parallel over draws, backend picked automatically for the platform
#' result <- big_b_and_M(h = 4, n_draws = 1000, n_var = 3, n_p = 2,
#'                       matrices = matrices, n_cores = 2)
#' big_b <- result[[1]]
#' big_M <- result[[2]]
#' }
#' @export
#' @import dplyr
big_b_and_M<-function(h,n_draws,n_var,n_p,data_=NULL,matrices=NULL,
                      n_cores=1,parallel=c("auto","fork","psock","none")){
  parallel <- match.arg(parallel)
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

  # ship only what the kernel needs (keeps PSOCK serialization light)
  mats <- list(B_list = matrices$B_list, M = matrices$M, intercept = matrices$intercept)

  per_draw <- apply_over_draws(n_draws = n_draws, n_cores = n_cores, parallel = parallel,
                               h = h, n_var = n_var, n_p = n_p, data_ = data_, mats = mats)

  big_b <- abind::abind(
    lapply(per_draw, function(x) array(x$b_h, dim = c(1, n_var * h, 1))),
    along = 3
  )

  big_M <- array(0, dim = c(n_var * h, n_var * h, n_draws))
  for (d in seq_len(n_draws)) {
    M_h_draw <- per_draw[[d]]$M_h
    for (cnt in 1:h) {
      zz <- 1
      for (cnt2 in cnt:h) {
        big_M[(1 + n_var * (cnt - 1)):(cnt * n_var),
              (1 + n_var * (cnt2 - 1)):(cnt2 * n_var),
              d] <- M_h_draw[[zz]]
        zz <- zz + 1
      }
    }
  }
  return(list(big_b,big_M))
}
