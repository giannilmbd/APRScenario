# setup.R
library(APRScenario)
library(zoo)
library(dplyr)



system("mkdir -p figures")
data("NKdata", package = "APRScenario")

# Rename the first column and process the data
names(NKdata)[1] <- "year"
NKdata$year <- zoo::as.yearqtr(NKdata$year, "%d-%b-%Y")
last_obs <- NKdata$year %>% dplyr::last()
varbls <<- names(NKdata)[-1]
X0 <<- NKdata[, varbls]

sr <<- matrix(c(
  1, -1, -1, # GDP
  1, 1, -1,  # inflation
  1, NA, 1   # interest rate
), nrow = 3, byrow = TRUE)

n_draws <<- 1500 # increase for final estimation
p <<- 4          # lags

# Possible subsampling
X <<- X0[1:(nrow(X0) - 4), ]

# Specify the model
specification <<- bsvarSIGNs::specify_bsvarSIGN$new(data = as.matrix(X, 3, 3),
                                                    p = p,
                                                    sign_irf = sr[1:3, 1:3])

# Estimate the model
posterior <<- bsvars::estimate(specification, S = n_draws)

# Compute and plot impulse responses
irf <<- bsvars::compute_impulse_responses(posterior, horizon = 40)
