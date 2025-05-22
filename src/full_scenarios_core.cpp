// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List full_scenarios_core(const arma::cube& big_b, const arma::cube& big_M,
                         const IntegerVector& obs, const NumericMatrix& path,
                         const IntegerVector& shocks, int h, int n_var,
                         Nullable<arma::vec> g_ = R_NilValue,
                         Nullable<arma::mat> Sigma_g_ = R_NilValue) {

  int n_draws = big_b.n_slices;
  int k_0 = obs.size() * h;
  bool has_shocks = std::any_of(shocks.begin(), shocks.end(), [](int val) {
    return val != NA_INTEGER;
  });
  int k_s = has_shocks ? shocks.size() * h : 0;

  vec g = g_.isNotNull() ? as<vec>(g_) : vec(k_s, fill::zeros);
  mat Sigma_g = Sigma_g_.isNotNull() ? as<mat>(Sigma_g_) : eye(k_s, k_s);
  mat L_g = chol(Sigma_g, "lower");

  std::vector<mat> mu_eps_vec(n_draws);
  std::vector<mat> Sigma_eps_vec(n_draws);
  std::vector<mat> mu_y_vec(n_draws);
  std::vector<mat> Sigma_y_vec(n_draws);

  for (int d = 0; d < n_draws; ++d) {
    mat b = big_b.slice(d);
    mat M = big_M.slice(d);
    mat M_t = trans(M);
    mat b_t = trans(b);

    mat C_h(k_0, n_var * h, fill::zeros);
    for (int j = 0; j < obs.size(); ++j) {
      int var_idx = obs[j] - 1;
      for (int t = 0; t < h; ++t) {
        C_h(j * h + t, var_idx + t * n_var) = 1.0;
      }
    }

    mat f(k_0 + k_s, 1, fill::zeros);
    for (int j = 0; j < obs.size(); ++j) {
      for (int t = 0; t < h; ++t) {
        int i_f = j * h + t;
        f(i_f, 0) = path(j, t);
      }
    }

    mat C_hat = C_h;
    if (has_shocks) {
      mat Xi(k_s, n_var * h, fill::zeros);
      for (int j = 0; j < shocks.size(); ++j) {
        if (shocks[j] == NA_INTEGER) continue;
        int var_idx = shocks[j] - 1;
        for (int t = 0; t < h; ++t) {
          Xi(j * h + t, var_idx + t * n_var) = 1.0;
        }
      }
      mat C_l = Xi * inv(M_t);
      C_hat = join_vert(C_h, C_l);
      f.rows(k_0, k_0 + k_s - 1) = C_l * b_t + g;
    }

    mat D = C_hat * M_t;
    mat D_ast = pinv(D);

    mat Omega_f = C_h * M_t * trans(C_h * M_t);
    if (has_shocks) {
      mat Z0(k_0, k_s, fill::zeros);
      mat Z1(k_s, k_0, fill::zeros);
      mat Omega_ext = join_vert(
        join_horiz(Omega_f, Z0),
        join_horiz(Z1, Sigma_g)
      );
      Omega_f = Omega_ext;
    }

    mat mu_eps = D_ast * (f - C_hat * b_t);
    mat I = eye(D_ast.n_rows, D_ast.n_rows);
    mat Sigma_eps = D_ast * Omega_f * D_ast.t() + (I - D_ast * D) * (I - D_ast * D).t();

    vec eps_draw = mu_eps;
    if (has_shocks) {
      vec z = randn(k_s);
      vec delta = g + L_g * z;
      for (int j = 0; j < shocks.size(); ++j) {
        if (shocks[j] == NA_INTEGER) continue;
        int var_idx = shocks[j] - 1;
        for (int t = 0; t < h; ++t) {
          int i_eps = var_idx + t * n_var;
          int i_delta = j * h + t;
          eps_draw(i_eps) = delta(i_delta);
        }
      }
    }

    mat mu_y = b_t + M_t * mu_eps;
    mat Sigma_y = M_t * M + (M_t * D_ast) * (Omega_f - D * D.t()) * (D_ast.t() * M);

    mu_eps_vec[d] = mu_eps;
    Sigma_eps_vec[d] = Sigma_eps;
    mu_y_vec[d] = mu_y;
    Sigma_y_vec[d] = Sigma_y;
  }

  return List::create(
    Named("mu_eps") = mu_eps_vec,
    Named("Sigma_eps") = Sigma_eps_vec,
    Named("mu_y") = mu_y_vec,
    Named("Sigma_y") = Sigma_y_vec
  );
}
