// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Draw-level ingredients of the Waggoner-Zha / APR (2021) Algorithm 1 posterior.
//
// For draw j the unconditional forecast is y_f = b_j + M_j' eps, eps ~ N(0, I),
// so log p(y_f | theta_j, Y) = -nh/2 log(2 pi) - log|det M_j| - 0.5 ||W_j (y_f - b_j)||^2
// with W_j = (M_j')^{-1}.

// Precompute W_j and log|det M_j| for every draw.
// [[Rcpp::export]]
List wz_prepare_core(const arma::cube& big_M) {
  const uword nh = big_M.n_rows, N = big_M.n_slices;
  cube W(nh, nh, N);
  vec logdet(N);
  for (uword j = 0; j < N; ++j) {
    mat Mt = trans(big_M.slice(j));
    W.slice(j) = inv(Mt);
    double val, sign;
    log_det(val, sign, Mt);
    logdet(j) = val;
  }
  return List::create(Named("W") = W, Named("logdet") = logdet);
}

// log p(y | theta_j, Y) for all draws j, up to the common constant.
static vec loglik_all(const vec& y, const mat& B, const cube& W, const vec& logdet) {
  const uword N = W.n_slices;
  vec ll(N);
  for (uword j = 0; j < N; ++j) {
    vec z = W.slice(j) * (y - B.col(j));
    ll(j) = -logdet(j) - 0.5 * dot(z, z);
  }
  return ll;
}

// Algorithm 1 of Antolin-Diaz, Petrella and Rubio-Ramirez (2021) with step 1
// carried out on a fixed set of draws from p(theta | Y): given a path y_f,
// p(theta_j | Y, y_f) is proportional to p(y_f | theta_j, Y) on that set.
// Step 2 draws y_f ~ N(mu_y_i, Sigma_y_i), the scenario distribution of draw i.
//
// b       nh x N   unconditional forecast means
// W, logdet        from wz_prepare_core()
// mu_y    nh x N   scenario means;  F  nh x nh x N, F_i F_i' = Sigma_y_i
// Returns the Rao-Blackwellised marginal weights pi (mean over iterations of
// the step-1 weight vector), the step-1 ESS at every iteration, and the chain.
// [[Rcpp::export]]
List wz_sirg_core(const arma::mat& b, const arma::cube& W, const arma::vec& logdet,
                  const arma::mat& mu_y, const arma::cube& F,
                  int n_iter, int burn, int start) {
  const uword N = W.n_slices, nh = b.n_rows;
  vec pi(N, fill::zeros);
  NumericVector ess(n_iter);
  IntegerVector chain(n_iter);
  uword idx = static_cast<uword>(start - 1);

  for (int it = 0; it < burn + n_iter; ++it) {
    if (it % 100 == 0) Rcpp::checkUserInterrupt();
    vec z = as<vec>(Rcpp::rnorm(nh));
    vec y = mu_y.col(idx) + F.slice(idx) * z;

    vec ll = loglik_all(y, b, W, logdet);
    vec w = exp(ll - ll.max());
    w /= accu(w);

    double u = R::unif_rand(), cum = 0.0;
    uword next = N - 1;
    for (uword j = 0; j < N; ++j) { cum += w(j); if (u <= cum) { next = j; break; } }
    idx = next;

    if (it >= burn) {
      pi += w;
      ess[it - burn] = 1.0 / accu(square(w));
      chain[it - burn] = static_cast<int>(idx) + 1;
    }
  }
  pi /= static_cast<double>(n_iter);
  return List::create(Named("pi") = pi, Named("ess") = ess, Named("chain") = chain);
}
