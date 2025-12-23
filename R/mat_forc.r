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
#' @param h (integer) forecast horison
#' @param n_draws (integer) Number of draws
#' @param n_var (integer) Number of variables
#' @param n_p (integer) Number of lags
#' @param data_ (matrix optional) The data, stacking Y over X (data and laggs)
#'        -- columns are observations (default taken from matrices$Z)
#'        NB: this is not necessarily the same as the data used to estimate the model
#'        If run counterfactuals in previoius historical period (ie not forecast) must pass the data up to previous period relative to counterfactual
#' @param matrices Optional matrices object from gen_mats() (default taken from calling environment)
#' @param max_cores maximum number of cores to use for parallel processing (default: 1 for Windows compatibility)
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
  # Default to 1 core for Windows compatibility
  cores <- 1
  
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
  # Derivation of the b and M matrices to produce the conditional MA form y=b+sum(epsilon*M)
  # See appendix of Antolin-Diaz et al (JME 21)
  ### b block
  # initial value is identity
  K_0<-array(rep(diag(n_var),n_draws),dim=c(n_var,n_var,n_draws))

  # NOTE that the facto K_{h-1} is in position K_[[h]] since K_{h} never appears in the algebra
  K_h=list()
  K_h[[1]]<-K_0

  if(h>1){
    for (i in 2:h){
      tmp2<-0
      for(j in 1:(i-1)){

        tmp2<-tmp2+parallel::mclapply(1:n_draws, function(d){
          if((j)<=n_p){
            tmp1<-matrices$B_list[[j]][,,d]
          }else{
            tmp1<-diag(0,n_var)
          }
          (K_h[[i-j]][,,d])%*%(tmp1)} , mc.cores = cores) %>% abind::abind(.,along=3)

      }
      K_h[[i]]<-tmp2+K_0

    }
  }


  M_h=list()
  M_h[[1]]<-matrices$M
  if(h>1){
    for(i in 2:h){
      tmp2<-0
      for(j in 1:min(i-1,n_p)){

        tmp2<-tmp2+parallel::mclapply(1:n_draws, function(d){M_h[[j]][,,d]%*%matrices$B_list[[j]][,,d]} , mc.cores = cores) %>% abind::abind(.,along=3)

      }
      M_h[[i]]<-tmp2
    }
  }


  # for each lag, start from the lag-matrix B then augment
  N_p_list<-list() # this is going to be a list of lists: for each lag and horizon

  for(l in 1:n_p){ # loop through each of the lags
    tmp00<-list()
    tmp00[[1]]<-matrices$B_list[[l]]
    if(h>1){ # if h>1 then need to update each of the N_s

      for (i in 2:h){
        tmp2<-0

        for(j in 1:min((i-1),n_p)){# new value is old value plus N(old value)*B

          tmp2<-tmp2+parallel::mclapply(1:n_draws, function(d){

            tmp00[[i-j]][,,d]%*%matrices$B_list[[j]][,,d]}
            ,mc.cores = cores
          ) %>% abind::abind(.,along=3)


        }
        if((l+i-1)<=n_p){ # add the linear term if exist
          tmp=matrices$B_list[[l+i-1]]
        }else{
          tmp=0
        }

        tmp00[[i]]<-tmp2+tmp
      }

    }
    N_p_list[[l]]<-tmp00 # lag is the first index of the list; horizon is the second
  }


  # b_h in two parts: the intercept and the lagged variables (from t-p to t)
  b_h<-parallel::mclapply(1:n_draws, function(d){matrices$intercept[,d]%*%K_h[[h]][,,d]} , mc.cores = cores) %>% abind::abind(.,along=3)

  ## part including lagged variables
  ### first create zeros 1xn_varxn_draws
  ylagged<-rep(0,n_var);
  ylagged<-array(rep(ylagged,n_draws),dim=c(1,n_var,n_draws))
  # lagged data

  for(cnt in 1:n_p){
    ylagged<-parallel::mclapply(1:n_draws,function(d){
      ylagged[,,d]+t(data_[(1+n_var*(cnt-1)):(n_var*cnt),ncol(data_)])%*%N_p_list[[cnt]][[h]][,,d]}, mc.cores = cores) %>% abind::abind(.,along=3)

  }
  
  # # for debugging
  # for(cnt in 1:n_p){
  #   ylagged<-lapply(1:n_draws,function(d){
  #     ylagged[,,d]+t(data_[(1+n_var*(cnt-1)):(n_var*cnt),ncol(data_)])%*%N_p_list[[cnt]][[h]][,,d]}) %>% abind::abind(.,along=3)
  #   
  # }
  
  # fix the order of dimensions of ylagged
  b_h<-b_h+(ylagged %>% array(.,dim=c(1,n_var,n_draws)))

  return(list(b_h,M_h))
}


