#' simscen function
#'
#' @param mu_eps mean innovation
#' @param Sigma_eps variance innovation
#' @param mu_y mean forecast
#' @param Sigma_y variance forecast
#' @param big_b history forecast
#' @param big_M IRF (innovation loading)
#' @param n_sim number of simulations
#' @param h horizon
#'
#' @returns conditional forecast path and distribution
#'
#' @export
#'
#' @import dplyr

SimScen<-function(mu_eps,Sigma_eps,mu_y,Sigma_y,big_b,big_M,n_sim,h,varbls){
 n_var=length(varbls)
  #draw shocks
  epsilon <- parallel::mclapply(1:n_draws, FUN = function(i) {
    # Generate n_sim samples for each draw
    samples <- MASS::mvrnorm(n = n_sim, mu = mu_eps[,1, i], Sigma = Sigma_eps[,, i])
    return(samples)
  }, mc.cores = parallel::detectCores()-1) %>% simplify2array() %>% aperm(.,c(2,1,3)) %>% array(.,dim=c(dim(big_b)[2],n_sim*n_draws))

  #expand big_b along simulations
  big_b_sim<-parallel::mclapply(1:n_draws,FUN=function(d){
    array(rep(big_b[1,,d],n_sim),dim=c(dim(big_b)[2],n_sim))},mc.cores=parallel::detectCores()-1) %>% simplify2array()%>% array(.,dim=c(dim(big_b)[2],n_sim*n_draws))
  # expand big_M
  big_M_sim<-parallel::mclapply(1:n_draws,FUN=function(d){
    array(rep(big_M[,,d],n_sim),dim=c(dim(big_M)[1],dim(big_M)[2],n_sim))
  },mc.cores=parallel::detectCores()-1) %>% simplify2array()%>% array(.,dim=c(dim(big_M)[1],dim(big_M)[2],n_sim*n_draws))
  # b+M'epsilon

  y_h<-parallel::mclapply(1:(n_draws*n_sim),FUN=function(d){
    big_b_sim[,d]+t(big_M_sim[,,d])%*%epsilon[,d]
  },mc.cores = parallel::detectCores()-1) %>% simplify2array() %>% array(.,dim=c(dim(big_M)[1],n_draws*n_sim))
  rownames(y_h)<-paste(rep(varbls,h),sort(rep(1:h,n_var)),sep='.')
  return(y_h)
}
