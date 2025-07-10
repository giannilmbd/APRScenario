# Declare global variables to avoid R CMD check NOTEs
# These are variables used by dplyr, ggplot2, and gen_mats global assignments

utils::globalVariables(c(
  # Data matrices created by gen_mats (global assignments with <<-)
  "Z", "Y", "XX", "B", "B_list", "M", "M_inv", "M_list", "intercept",
  # Posterior object and specification
  "posterior", "specification", "n_var", "n_p",
  # dplyr/tidyr variables
  ".", "variable", "hor", "data", "density", "time", "vars", "h",
  # ggplot2 aesthetic variables  
  "lower_bound", "upper_bound", "KLM", "center", "lower", "upper", "hist",
  "Median", "LB", "UB", "LB_s", "UB_s"
))
