% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scenarios_old2.r
\name{scenarios_old2}
\alias{scenarios_old2}
\title{scenarios function (optimized with Rcpp)}
\usage{
scenarios_old2(
  h = 3,
  path = NULL,
  obs = NULL,
  shocks = NULL,
  n_draws,
  n_sample = n_draws,
  n_var,
  n_p,
  data_ = Z
)
}
\arguments{
\item{h}{forecast horizon}

\item{path}{conditional path of observables}

\item{obs}{position of observable(s)}

\item{shocks}{position of non-driving shocks (NA if all driving)}

\item{n_draws}{Number of draws}

\item{n_sample}{Number of draws to sample (must be <= n_draws)}

\item{n_var}{Number of variables}

\item{n_p}{Number of lags}

\item{data_}{optional matrix of Y over X (default Z)}
}
\value{
list of mu_eps, Sigma_eps, mu_y, Sigma_y, big_b, big_M, draws_used
}
\description{
scenarios function (optimized with Rcpp)
}
