#' big_b_and_M
#' This function returns the extended b and M matrices as in APR
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
#' @returns the big_b and big_M matrices of mean and IRF
#' @examples
#' \dontrun{
#' # Example usage for creating extended matrices
#' result <- big_b_and_M(h = 4, n_draws = 1000, n_var = 3, n_p = 2,
#'                       matrices = matrices)
#' big_b <- result[[1]]
#' big_M <- result[[2]]
#' }
#' @export
#' @import dplyr
big_b_and_M<-function(h,n_draws,n_var,n_p,data_=NULL,matrices=NULL){
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
  big_b<-array(0,dim=c(1,n_var*h,n_draws))
  big_M<-array(0,dim=c(n_var*h,n_var*h,n_draws))
  M_h=list()
  # Initialize M_h with the first M matrix
  for(cnt in 1:h){
    tmp<-mat_forc(h = cnt, n_draws = n_draws, n_var = n_var, n_p = n_p, data_ = data_)
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
