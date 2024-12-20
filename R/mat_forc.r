#' mat_forc function
#'##############################################################################
#'  NB: HERE WE USE Antolin-Diaz et al notation                                #
#'  B is reduced form;                                                         #
#'  A is structural;                                                           #
#'  d is intercepts                                                            #
#'  M is reduced so that E(u*u')=Sigma=(A_0*A_0')^{-1} and M_0=A_0^{-1}*Q      #
#'  Note that the code returns conflicting notation:                           #
#'  B=>A_0^{-1}*Q and                                                          #
#'  A=>B                                                                       #
#'##############################################################################
#' @param h forecast horison
#'
#' @returns the big_b and big_M matrices of mean and IRF
#' @export
#' @import dplyr
#'
mat_forc<-function(h=1){
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
            tmp1<-B_list[[j]][,,d]
          }else{
            tmp1<-diag(0,n_var)
          }
          (K_h[[i-j]][,,d])%*%(tmp1)} , mc.cores = parallel::detectCores()-1) %>% simplify2array()

      }
      K_h[[i]]<-tmp2+K_0

    }
  }


  M_h=list()
  M_h[[1]]<-M
  if(h>1){
    for(i in 2:h){
      tmp2<-0
      for(j in 1:min(i-1,n_p)){

        tmp2<-tmp2+parallel::mclapply(1:n_draws, function(d){M_h[[j]][,,d]%*%B_list[[j]][,,d]} , mc.cores = parallel::detectCores()-1) %>% simplify2array()

      }
      M_h[[i]]<-tmp2
    }
  }


  # for each lag, start from the lag-matrix B then augment
  N_p_list<-list() # this is going to be a list of lists: for each lag and horizon

  for(l in 1:n_p){ # loop through each of the lags
    tmp00<-list()
    tmp00[[1]]<-B_list[[l]]
    if(h>1){ # if h>1 then need to update each of the N_s

      for (i in 2:h){
        tmp2<-0

        for(j in 1:min((i-1),n_p)){# new value is old value plus N(old value)*B

          tmp2<-tmp2+parallel::mclapply(1:n_draws, function(d){

            tmp00[[i-j]][,,d]%*%B_list[[j]][,,d]}
            ,mc.cores = parallel::detectCores()-1
          ) %>% simplify2array()


        }
        if((l+i-1)<=n_p){ # add the linear term if exist
          tmp=B_list[[l+i-1]]
        }else{
          tmp=0
        }

        tmp00[[i]]<-tmp2+tmp
      }

    }
    N_p_list[[l]]<-tmp00 # lag is the first index of the list; horizon is the second
  }


  # b_h in two parts: the intercept and the lagged variables (from t-p to t)
  b_h<-parallel::mclapply(1:n_draws, function(d){intercept[,d]%*%K_h[[h]][,,d]} , mc.cores = parallel::detectCores()-1) %>% simplify2array(.,higher=T)

  ## part including lagged variables
  ### first create zeros 1xn_varxn_draws
  ylagged<-rep(0,n_var);
  ylagged<-array(rep(ylagged,n_draws),dim=c(1,n_var,n_draws))
  # lagged data

  for(cnt in 1:n_p){
    ylagged<-parallel::mclapply(1:n_draws,function(d){
      ylagged[,,d]+t(Z[(1+n_var*(cnt-1)):(n_var*cnt),ncol(Z)])%*%N_p_list[[cnt]][[h]][,,d]}, mc.cores = parallel::detectCores()-1) %>% simplify2array()

  }
  # fix the order of dimensions of ylagged
  b_h<-b_h+(ylagged %>% array(.,dim=c(1,n_var,n_draws)))

  return(list(b_h,M_h))
}

big_b_and_M<-function(h){
  big_b<-array(0,dim=c(1,n_var*h,n_draws))
  big_M<-array(0,dim=c(n_var*h,n_var*h,n_draws))
  M_h=list()
  for(cnt in 1:h){
    tmp<-mat_forc(cnt)
    b_h<-tmp[[1]]
    M_h[[cnt]]<-tmp[[2]][[cnt]]
    big_b[1,(1+n_var*(cnt-1)):(cnt*n_var),]<-b_h

  }
  for(cnt in 1:h){
    zz=1
    for(cnt2 in cnt:h){
      big_M[(1+n_var*(cnt-1)):(cnt*n_var),(1+n_var*(cnt2-1)):(cnt2*n_var),]<-M_h[[zz]]
      zz=zz+1
    }
  }
  return(list(big_b,big_M))
}
