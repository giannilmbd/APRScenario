% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forc_h.r
\name{forc_h}
\alias{forc_h}
\title{forc_h function}
\usage{
forc_h(h = 1, n_sim = 200, n_draws, n_var, n_p, data_ = Z)
}
\arguments{
\item{h}{forecast horizon}

\item{n_sim}{length of shock simulation}

\item{n_draws}{Number of draws}

\item{n_var}{Number of variables}

\item{n_p}{Number of lags}

\item{data_}{Optional matrix of data n_var*h+1 x T (default is Z)}
}
\value{
a matrix of unconditional forecasts
}
\description{
forc_h function
}
