// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <iostream>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List scenarios_core(List big_b_list, List big_M_list, List C_hat_list, List D_list, List D_ast_list, 
                    List f_hat_list, List Omega_hat_list) {
  int n_draws = big_b_list.size();

  std::vector<arma::mat> mu_eps_vec(n_draws);
  std::vector<arma::mat> Sigma_eps_vec(n_draws);
  std::vector<arma::mat> mu_y_vec(n_draws);
  std::vector<arma::mat> Sigma_y_vec(n_draws);
  vec      progress = round(linspace(0, n_draws, 10));
  Rcpp::Rcout << "Starting scenario computations...\n" << std::flush;
  Progress p(10, true);  // true = show progress bar in R
  

  for (int d = 0; d < n_draws; ++d) {
    
    if (Progress::check_abort()) return List::create();  // allow interrupt
    
    arma::mat big_b  = big_b_list[d];
    arma::mat big_M  = big_M_list[d];
    arma::mat C_hat  = C_hat_list[d];
    arma::mat D      = D_list[d];
    arma::mat D_ast  = D_ast_list[d];
    arma::mat f_hat  = f_hat_list[d];
    arma::mat Omega_hat = Omega_hat_list[d];

    // mu_eps = D_ast %*% (f_hat - C_hat %*% big_b)
    arma::mat mu_eps = D_ast * (f_hat - C_hat * big_b);

    // Sigma_eps
    arma::mat I = eye(D_ast.n_rows, D_ast.n_rows);
    arma::mat Sigma_eps = D_ast * Omega_hat * D_ast.t() +
                          (I - D_ast * D) * (I - D_ast * D).t();

    // mu_y = big_b + t(big_M) %*% mu_eps
    arma::mat mu_y = big_b + big_M.t() * mu_eps;

    // Sigma_y
    arma::mat Sigma_y = big_M.t() * big_M +
                        (big_M.t() * D_ast) * (Omega_hat - D * D.t()) * (D_ast.t() * big_M);

    mu_eps_vec[d] = mu_eps;
    Sigma_eps_vec[d] = Sigma_eps;
    mu_y_vec[d] = mu_y;
    Sigma_y_vec[d] = Sigma_y;

    if (any(progress == d)) {
      p.increment();
    }
  }

  return List::create(
    Named("mu_eps") = mu_eps_vec,
    Named("Sigma_eps") = Sigma_eps_vec,
    Named("mu_y") = mu_y_vec,
    Named("Sigma_y") = Sigma_y_vec
  );
}
// install.packages("Rcpp")
// install.packages("RcppArmadillo")
// library(Rcpp)
// library(RcppArmadillo)
// Rcpp::sourceCpp("scenarios_core.cpp")
