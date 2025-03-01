% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_cond_histo.r
\name{plot_cond_histo}
\alias{plot_cond_histo}
\title{plot_cond_histo function}
\usage{
plot_cond_histo(variable = NULL, horizon = 1, threshold = NULL, data = NULL)
}
\arguments{
\item{variable}{Name of variable to be plotted (string)}

\item{horizon}{At which horizon (horizon<=h)}

\item{threshold}{(optional) If present compute P(x>threshold)}

\item{data}{data of conditional forecasts}
}
\value{
ggplot object (plot)
}
\description{
This function uses the conditional probability calculations (eg scenarios) and plots the histogram of the selected variable
}
