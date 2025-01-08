#' KL function
#' APR suggest this measure to assess the "plausibility" of the conditional forecast.
#' It is based on the Kullback-Leibler measure of distance between the unconditional forecast and the
#' conditiona/scenario forecast.
#'
#' @param mu_eps mean of innovation
#' @param Sigma_eps variance of innovation
#' @param plot_ logical; if TRUE then a histogram of the KL measure is returned
#'
#' @returns Returns the APR 'q': ie distance from a fair binomial distribution
#'
#' @export
#' @import dplyr
KL<-function(Sigma_eps,mu_eps,h,plot_=F){
  DKL<-parallel::mclapply(1:n_draws,FUN=function(d){
    0.5*(psych::tr(Sigma_eps[,,d])+t(mu_eps[,,d])%*%mu_eps[,,d]-n_var*h-log(det(Sigma_eps[,,d])))
  },mc.cores = parallel::detectCores()-1) %>% simplify2array()
  q<-0.5*(1+sqrt(1-exp(-2*DKL/h*n_var)))
  p=NA
  if(plot_){
    tmp<-data.frame(drw=seq_along(q),KLM=q)
    p<-ggplot(data=tmp)+geom_histogram(aes(x=KLM),alpha=.5)+
      geom_vline(xintercept = 0.5)+xlab('')+ylab('')+labs(title='Kullback-Leibler plausibility measure')+
      theme_minimal()
  }
  return(list(q,p))
}
