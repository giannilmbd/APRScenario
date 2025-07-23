#' gen_mats function
#'
#' this function returns the matrices necessary for forecasts
#'
#' @param posterior Posterior estimation results (eg from BsvarSIGNs)
#' @param specification Optional specification object (default taken from calling environment)
#' @returns Returns all objects necessary for scenario analysis (e.g., IRF matrix),
#'   including: \code{M}, \code{M_inv}, \code{M_list}, \code{B}, \code{B_list},
#'   \code{n_p}, \code{n_var}, \code{Y}, \code{X}, and \code{Z}.
#' @examples
#' \donttest{
#' library(APRScenario)
#' data(NKdata)
#'
#' # Minimal example with a toy specification
#' spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[,2:4]), p = 1)
#' est <- bsvars::estimate(spec, S = 10)  # Use small S for fast test
#' gen_mats(posterior = est, specification = spec)
#' }

#' @export
#'
#' @import dplyr



gen_mats<-function(posterior=NULL, specification=NULL){
  # Get specification from calling environment if not provided
  if(is.null(specification)) {
    if(exists("specification", envir=parent.frame())) {
      specification <- get("specification", envir=parent.frame())
    } else {
      stop("Please provide specification object or ensure it exists in calling environment")
    }
  }
  n_var<-dim(posterior$posterior$B)[1]
  n_p<-(dim(posterior$posterior$A)[2]-1)/n_var
  n_draws<-dim(posterior$posterior$B)[3]
  # this function returns the matrices necessary for forecasts
  ##### FIRST DEFINE NECESSARY MATRICES AND DIMENSIONS
  ## Extract matrices
  M_inv<-posterior$posterior$B # This already include Q i.e. B^{-1}=A_0^{-1}%*%Q
  #Q<-posterior$posterior$Q # Rotation matrices
  # remember to transpose
  # Use single core if running in CRAN environment or on Windows
  cores <- if(!identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "")) 1 else min(2, parallel::detectCores()-1)
  if (.Platform$OS.type == "windows") cores <- 1
  
  M_list <- parallel::mclapply(1:dim(M_inv)[3], function(d) solve(t(M_inv[,,d])), mc.cores = cores)
  M <-  abind::abind(M_list, along = 3)
  B<-posterior$posterior$A # this is lags plus constant (in the last row... see last row of specification$data_matrices$X)
  intercept<-B[,dim(B)[2],]
  B<-B[,-dim(B)[2],]

  # define dimensions
  n_var<-dim(B)[1]
  n_p<-(dim(B)[2])/n_var
  n_draws<-dim(B)[3]

  ## Extract lag-matrix-coefficients: from n_var \times n_var*p+1 \times n_draws to  n_var \times n_draws separate n_var*p matrices
  # REMEMBER TO TRANSPOSE
  B_list<-list()
  for(cnt in 1:n_p){
    tmp<-parallel::mclapply(1:n_draws,FUN = function(d){t(B[,(n_var*(cnt-1)+1):(n_var*cnt),d])},
                            mc.cores=cores) %>% simplify2array()
    B_list[[cnt]]<-tmp
  }
  # intercept parameters


  Y<-specification$data_matrices$Y
  XX<-specification$data_matrices$X

  Z<-rbind(Y,XX[-nrow(XX),])

  ######################################

  # Return named list with all necessary objects
  return(list(
    M = M,
    M_inv = M_inv,
    M_list = M_list,
    B = B,
    B_list = B_list,
    intercept = intercept,
    n_var = n_var,
    n_p = n_p,
    Y = Y,
    X = XX,
    Z = Z
  ))
}
