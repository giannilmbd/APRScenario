#' KL function
#' APR suggest this measure to assess the "plausibility" of the conditional forecast.
#' It is based on the Kullback-Leibler measure of distance between the unconditional forecast and the
#' conditional/scenario forecast.
#'
#' @param Sigma_eps variance of innovation
#' @param mu_eps mean of innovation
#' @param h forecast horizon
#' @param plot_ logical; if TRUE then a histogram of the KL measure is returned
#' @param max_cores maximum number of cores to use for parallel processing (default: NULL, uses CRAN-compliant detection with Windows=1)
#'
#' @returns Returns the APR 'q': ie distance from a fair binomial distribution
#' @examples
#' # Example with simulated innovation data
#' # Set dimensions
#' n_var <- 3
#' h <- 4
#' n_draws <- 10
#' n_innovations <- n_var * h
#' 
#' # Create simulated innovation means and covariances
#' set.seed(123)
#' mu_eps <- array(rnorm(n_innovations * 1 * n_draws, mean = 0, sd = 0.1), 
#'                 dim = c(n_innovations, 1, n_draws))
#' 
#' Sigma_eps <- array(0, dim = c(n_innovations, n_innovations, n_draws))
#' for (d in 1:n_draws) {
#'   temp_cov <- matrix(rnorm(n_innovations^2), n_innovations, n_innovations)
#'   Sigma_eps[,,d] <- temp_cov %*% t(temp_cov) + diag(n_innovations) * 0.5
#' }
#' 
#' # Calculate KL measure
#' kl_result <- KL(Sigma_eps, mu_eps, h, plot_ = FALSE)
#' print(head(kl_result[[1]]))  # Print first few q values
#'
#' @export
#' @import dplyr
#' @importFrom ggplot2 ggplot geom_histogram aes geom_vline xlab ylab labs theme_minimal
KL<-function(Sigma_eps,mu_eps,h,plot_=FALSE,max_cores=NULL){
  # Default to 1 core for Windows compatibility
  if(is.null(max_cores)) {
    cores <- 1
  } else {
    cores <- max_cores
  }
  
  n_var<-dim(Sigma_eps)[1]/h
  n_draws<-dim(Sigma_eps)[3]
  DKL<-parallel::mclapply(1:n_draws,FUN=function(d){
    0.5*(psych::tr(Sigma_eps[,,d])+t(mu_eps[,,d])%*%mu_eps[,,d]-n_var*h-log(det(Sigma_eps[,,d])))
  },mc.cores = cores) %>% simplify2array()

  q<-0.5*(1+sqrt(1-exp(-2*DKL/(h*n_var))))

  p=NA
  if(plot_){
    tmp<-data.frame(drw=seq_along(q),KLM=q)
    p<-ggplot(data=tmp)+geom_histogram(aes(x=KLM),alpha=.5)+
      geom_vline(xintercept = 0.5)+xlab('')+ylab('')+labs(title='Kullback-Leibler plausibility measure')+
      theme_minimal()
  }
  return(list(q,p))
}
