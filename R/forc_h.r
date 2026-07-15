#' forc_h function
#'
#' @param h forecast horizon
#' @param n_sim length of shock simulation
#' @param data_ Optional matrix of data n_var*h+1 x T. If NULL, defaults to matrices$Z
#' @param posterior Optional posterior object (default taken from calling environment)
#' @param matrices Optional matrices object from gen_mats() (default taken from calling environment)
#' @param max_cores number of workers used to parallelize over posterior draws
#'        (default 1 = serial; forked processes on Unix/macOS, PSOCK cluster on Windows;
#'        the shock-simulation loop additionally uses forking where available)
#' @returns a matrix of unconditional forecasts
#' @examples
#' \dontrun{
#' # Example usage for unconditional forecasting
#' forecast <- forc_h(h = 4, n_sim = 1000, 
#'                    posterior = posterior, matrices = matrices)
#' }
#' @export
#' @import dplyr

forc_h<-function(h=1,n_sim=200,data_=NULL,posterior=NULL,matrices=NULL,max_cores=1){
  # Get matrices from calling environment if not provided
  if(is.null(matrices)) {
    if(exists("matrices", envir=parent.frame())) {
      matrices <- get("matrices", envir=parent.frame())
    } else {
      stop("Please provide matrices object from gen_mats() or ensure it exists in calling environment")
    }
  }
  
  # Get posterior from calling environment if not provided
  if(is.null(posterior)) {
    if(exists("posterior", envir=parent.frame())) {
      posterior <- get("posterior", envir=parent.frame())
    } else {
      stop("Please provide posterior object or ensure it exists in calling environment")
    }
  }
  
  # Use matrices$Z as default for data_ only if not explicitly provided
  if(is.null(data_)) {
    data_ <- matrices$Z
  }
  
  n_var<-dim(posterior$posterior$B)[1]
  n_p<-(dim(posterior$posterior$A)[2]-1)/n_var
  n_draws<-dim(posterior$posterior$B)[3]

  # this function takes the forecast matrices b and M and combines with simulated future shocks
  # the output will be the forecsast on n_var, h-horizon,n_draws and n_sim (last two dims eventually mixed together)

  y_h=array(rep(rep(rep(rep(0,n_var),h),n_draws),n_sim),dim=c(n_var,h,n_draws,n_sim))
  # record also the shock and history parts separately
  shock_h=array(rep(rep(rep(rep(0,n_var),h),n_draws),n_sim),dim=c(n_var,h,n_draws,n_sim))

  hist_h=array(rep(rep(rep(rep(0,n_var),h),n_draws),n_sim),dim=c(n_var,h,n_draws,n_sim))

  # mclapply forks; on Windows only 1 core is supported
  cores <- if (.Platform$OS.type == "windows") 1L else max(1L, as.integer(max_cores))

  epsilon<-parallel::mclapply(1:n_draws,function(d)MASS::mvrnorm(n = n_sim, mu = rep(0,n_var), Sigma = diag(1,n_var)),
                              mc.cores = cores) %>% simplify2array()

  epsilon<-parallel::mclapply(1:h,function(i)epsilon,mc.cores=cores) %>% simplify2array()  %>%
    aperm(.,c(2,4,3,1))

  # forecast-mean terms and IRFs for all horizons in one pass per draw
  mats <- list(B_list = matrices$B_list, M = matrices$M, intercept = matrices$intercept)
  per_draw <- apply_over_draws(n_draws = n_draws, n_cores = max_cores, parallel = "auto",
                               h = h, n_var = n_var, n_p = n_p, data_ = data_, mats = mats)
  b_all <- vapply(per_draw, function(x) matrix(x$b_h, n_var, h), matrix(0, n_var, h)) # n_var x h x n_draws
  M_h <- lapply(seq_len(h), function(i)
    abind::abind(lapply(per_draw, function(x) x$M_h[[i]]), along = 3))

  for(cnt in 1:h){
    b_h<-array(b_all[, cnt, ], dim = c(1, n_var, n_draws))

    fut_shocks<-array(rep(rep(rep(0,n_var),n_draws),n_sim),dim=c(1,n_var,n_draws,n_sim))
    # for each step in h have to add up all future shocks
    for (tt in 1:cnt){
      tmp <- parallel::mclapply(1:n_sim, FUN = function(s) {
        # For each simulation, apply across draws
        sapply(1:n_draws, function(d) {
          # shock at future period tt hits y at period cnt via the (cnt-tt)-step IRF
          as.vector(epsilon[, tt, d, s] %*% M_h[[cnt - tt + 1]][,, d])
        })
      }, mc.cores = cores)

      # Convert the result into an array
      tmp <- simplify2array(tmp)

      fut_shocks<-fut_shocks+ array(tmp, dim = c(1, n_var, n_draws, n_sim))
    }
    shock_h[,cnt,,]<-fut_shocks[1,,,]

    hist_h[,cnt,,]<-array(rep(b_h,n_sim),dim=c(n_var,n_draws,n_sim))
    # add up b and M part (shocks)

    y_h[,cnt,,]<-array(rep(b_h,n_sim),dim=c(n_var,n_draws,n_sim))+ fut_shocks[1,,,]

  }
  shock_h_flattened<-array(shock_h,dim=c(n_var,h,n_draws*n_sim))
  hist_h_flattened<-array(hist_h,dim=c(n_var,h,n_draws*n_sim))
  y_h_flattened <- array(y_h, dim = c(n_var,h, n_draws * n_sim))
  return(list(y_h_flattened,shock_h_flattened,hist_h_flattened,b_h))
}
