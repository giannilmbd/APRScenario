#' gen_mats function
#'
#' this function assign to the global environment the matrices necessary for forecasts
#'
#' @returns Returns all the objects necessary for the other functions in scenario analysis (eg IRF matrix)
#' @export
#'
#' @import dplyr



gen_mats<-function(){
  # this function assign to the global environment the matrices necessary for forecasts
  ##### FIRST DEFINE NECESSARY MATRICES AND DIMENSIONS
  ## Extract matrices
  M_inv<<-posterior$posterior$B # This already include Q i.e. B^{-1}=A_0^{-1}%*%Q
  #Q<<-posterior$posterior$Q # Rotation matrices
  # remember to transpose
  M_list <<- parallel::mclapply(1:dim(M_inv)[3], function(d) solve(t(M_inv[,,d])), mc.cores = parallel::detectCores()-1)
  M <<- simplify2array(M_list)
  B<<-posterior$posterior$A # this is lags plus constant (in the last row... see last row of specification$data_matrices$X)
  intercept<<-B[,dim(B)[2],]
  B<<-B[,-dim(B)[2],]

  # define dimensions
  n_var<<-dim(B)[1]
  n_p<<-(dim(B)[2])/n_var
  # n_draws=dim(B)[3]

  ## Extract lag-matrix-coefficients: from n_var \times n_var*p+1 \times n_draws to  n_var \times n_draws separate n_var*p matrices
  # REMEMBER TO TRANSPOSE
  B_list<<-list()
  for(cnt in 1:n_p){
    tmp<-parallel::mclapply(1:n_draws,FUN = function(d){t(B[,(n_var*(cnt-1)+1):(n_var*cnt),d])},
                            mc.cores=parallel::detectCores()-1) %>% simplify2array()
    B_list[[cnt]]<<-tmp
  }
  # intercept parameters


  Y<<-specification$data_matrices$Y
  XX<<-specification$data_matrices$X

  Z<<-rbind(Y,XX[-nrow(XX),])

  ######################################

  # assign("M_list", M_list, envir = .GlobalEnv)
  # assign("B_list", B_list, envir = .GlobalEnv)
  # assign("Y", Y, envir = .GlobalEnv)
  # assign("Z", Z, envir = .GlobalEnv)
  # assign("n_var", n_var, envir = .GlobalEnv)
  # assign("n_p", n_p, envir = .GlobalEnv)
  # assign("M", M, envir = .GlobalEnv)
  # assign("intercept",intercept,envir=.GlobalEnv)
}
