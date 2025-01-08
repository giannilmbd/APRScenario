#' forc_h function
#'
#' @param h forecast horizon
#' @param n_sim length of shock simulation
#' @returns a matrix of unconditional forecasts
#' @export
#' @import dplyr

forc_h<-function(h=1,n_sim=200){

  # this function takes the forecast matrices b and M and combines with simulated future shocks
  # the output will be the forecsast on n_var, h-horizon,n_draws and n_sim (last two dims eventually mixed together)

  y_h=array(rep(rep(rep(rep(0,n_var),h),n_draws),n_sim),dim=c(n_var,h,n_draws,n_sim))
  # record also the shock and history parts separately
  shock_h=array(rep(rep(rep(rep(0,n_var),h),n_draws),n_sim),dim=c(n_var,h,n_draws,n_sim))
  hist_h=array(rep(rep(rep(rep(0,n_var),h),n_draws),n_sim),dim=c(n_var,h,n_draws,n_sim))

  epsilon<-parallel::mclapply(1:n_draws,function(d)MASS::mvrnorm(n = n_sim, mu = rep(0,n_var), Sigma = diag(1,n_var)),
                              mc.cores = parallel::detectCores()-1) %>% simplify2array()

  epsilon<-parallel::mclapply(1:h,function(i)epsilon,mc.cores=parallel::detectCores()-1) %>% simplify2array()  %>%
    aperm(.,c(2,4,3,1))


  for(cnt in 1:h){
    tmp<-mat_forc(cnt)
    b_h<-tmp[[1]]
    M_h<-tmp[[2]]

    fut_shocks<-array(rep(rep(rep(0,n_var),n_draws),n_sim),dim=c(1,n_var,n_draws,n_sim))
    # for each step in h have to add up all future shocks
    for (tt in 1:cnt){
      tmp <- parallel::mclapply(1:n_sim, FUN = function(s) {
        # For each simulation, apply across draws
        sapply(1:n_draws, function(d) {
          # Each draw returns a n_var vector
          as.vector(epsilon[, tt, d, s] %*% M_h[[tt]][,, d])
        })
      }, mc.cores = parallel::detectCores()-1)

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
