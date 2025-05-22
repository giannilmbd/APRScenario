// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List full_scenarios_core(const arma::cube& big_b, const arma::cube& big_M,
                         const IntegerVector& obs, const NumericVector& path,
                         const IntegerVector& shocks, int h, int n_var) {

  int n_draws = big_b.n_slices;
  int k_0 = obs.size() * h;

  // check if *any* value in shocks is not NA
  bool has_shocks = std::any_of(shocks.begin(), shocks.end(), [](int val) {
    return val != NA_INTEGER;
  });

  int k_s = has_shocks ? shocks.size() * h : 0;

  // Input validation
  if (path.size() < k_0) {
    stop("Path vector is too short. Expected length %d, got %d", k_0, path.size());
  }

  std::vector<mat> mu_eps_vec(n_draws);
  std::vector<mat> Sigma_eps_vec(n_draws);
  std::vector<mat> mu_y_vec(n_draws);
  std::vector<mat> Sigma_y_vec(n_draws);

  // Progress p(n_draws, true);
  Rcpp::Rcout << "Running scenarios...\n" << std::flush;

  for (int d = 0; d < n_draws; ++d) {
    // if (Progress::check_abort()) return List::create();

    mat b = big_b.slice(d);   // (1 x k)
    mat M = big_M.slice(d);   // (nM x k)
    mat M_t = trans(M);
    mat b_t = trans(b);

    // Build C_h
    mat C_h(k_0, n_var * h, fill::zeros);
    for (int j = 0; j < obs.size(); ++j) {
      int var_idx = obs[j] - 1;
      for (int t = 0; t < h; ++t) {
        C_h(j * h + t, var_idx + t * n_var) = 1.0;
      }
    }

    mat f(k_0 + k_s, 1, fill::zeros);
    for (int i = 0; i < k_0; ++i) {
      f(i, 0) = path[i];
    }

    mat C_hat = C_h;

    if (has_shocks) {
      mat Xi(k_s, n_var * h, fill::zeros);
      for (int j = 0; j < shocks.size(); ++j) {
        if (shocks[j] == NA_INTEGER) continue;  // skip if still NA
        int var_idx = shocks[j] - 1;
        for (int t = 0; t < h; ++t) {
          Xi(j * h + t, var_idx + t * n_var) = 1.0;
        }
      }

      mat C_l = Xi * inv(M_t);
      C_hat = join_vert(C_h, C_l);

      f.rows(k_0, k_0 + k_s - 1) = C_l * b_t;
    }

    mat D = C_hat * M_t;
    mat D_ast = pinv(D);

    mat Omega_f = C_h * M_t * trans(C_h * M_t);
    mat Omega_hat = Omega_f;
    if (has_shocks) {
      mat Z0(k_0, k_s, fill::zeros);
      mat Z1(k_s, k_0, fill::zeros);
      mat Is = eye(k_s, k_s);
      Omega_hat = join_vert(
        join_horiz(Omega_f, Z0),
        join_horiz(Z1, Is)
      );
    }

    mat mu_eps = D_ast * (f - C_hat * b_t);

    mat I = eye(D_ast.n_rows, D_ast.n_rows);
    mat Sigma_eps = D_ast * Omega_hat * D_ast.t() +
                    (I - D_ast * D) * (I - D_ast * D).t();

    mat mu_y = b_t + M_t * mu_eps;
    mat Sigma_y = M_t * M +
                  (M_t * D_ast) * (Omega_hat - D * D.t()) * (D_ast.t() * M);

    mu_eps_vec[d] = mu_eps;
    Sigma_eps_vec[d] = Sigma_eps;
    mu_y_vec[d] = mu_y;
    Sigma_y_vec[d] = Sigma_y;

    // p.increment();
  }

  return List::create(
    Named("mu_eps") = mu_eps_vec,
    Named("Sigma_eps") = Sigma_eps_vec,
    Named("mu_y") = mu_y_vec,
    Named("Sigma_y") = Sigma_y_vec
  );
}
