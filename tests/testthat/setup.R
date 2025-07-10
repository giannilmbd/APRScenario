# setup.R
library(APRScenario)
library(zoo)
library(dplyr)

# Load test data
data("NKdata", package = "APRScenario")

# Basic setup for testing - no model estimation
varbls <<- names(NKdata)[-1]
X0 <<- NKdata[, varbls]

# Create minimal test matrices for package functions
n_var <<- 3
n_p <<- 4
n_draws <<- 10  # minimal for testing
