#' plot_bvars: This function plots the IRFs generated with the BVAR
#'
#'
#'
#' @param M IRFs produced by eg bvarSIGNs
#' @param significance_level (eg 0.05)
#' @param central_tendency eg 'mean' or 'median'
#' @param variable_names vector of names of variables (strings)
#' @param shock_names vector of names of variables (strings)
#'
#' @returns a list of ggplot objects (plots)
#' @examples
#' # Example with simulated IRF data
#' # Create simulated IRF array (n_vars, n_shocks, n_periods, n_draws)
#' set.seed(123)
#' n_vars <- 3
#' n_shocks <- 3
#' n_periods <- 10
#' n_draws <- 50
#' 
#' # Generate IRF responses that decay over time
#' M <- array(0, dim = c(n_vars, n_shocks, n_periods, n_draws))
#' for (i in 1:n_vars) {
#'   for (j in 1:n_shocks) {
#'     for (t in 1:n_periods) {
#'       # Create decaying responses with some randomness
#'       base_response <- ifelse(i == j, 1, 0.3) * exp(-0.1 * (t-1))
#'       M[i, j, t, ] <- rnorm(n_draws, mean = base_response, sd = 0.1)
#'     }
#'   }
#' }
#' 
#' # Create plots
#' var_names <- c("GDP", "CPI", "FFR")
#' shock_names <- c("Supply", "Demand", "Monetary")
#' 
#' plots <- plot_bvars(M, 
#'                     variable_names = var_names,
#'                     shock_names = shock_names,
#'                     significance_level = 0.1)
#' 
#' # plots is a list of ggplot objects
#' print(length(plots))
#'
#' @export
#' @import dplyr
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon xlab ylab labs theme element_text geom_hline element_blank element_line
#' @importFrom stats median quantile

plot_bvars <- function(M, significance_level = 0.05, central_tendency = 'mean', variable_names = NULL, shock_names = NULL) {


  # Get dimensions of the array
  N <- dim(M)[1]  # Number of variables
  s <- dim(M)[2]  # Number of shocks
  T <- dim(M)[3]  # Number of periods
  D <- dim(M)[4]  # Number of draws

  # Set default variable and shock names if not provided
  if (is.null(variable_names)) {
    variable_names <- paste('Variable', 1:N)
  }
  if (is.null(shock_names)) {
    shock_names <- paste('Shock', 1:s)
  }

  # Initialize a list to store plots
  plot_list <- list()
  plot_counter <- 1

  # Loop over each variable and shock
  for  (j in 1:s){
    for (i in 1:N) {
      # Extract the T x D matrix for the current variable and shock
      # M_ij is a T x D matrix where each row corresponds to a time period
      # and each column corresponds to a draw
      M_ij <- M[i, j, , ]  # T x D matrix

      # Compute the central tendency (mean or median) across draws for each time period
      if (central_tendency == 'mean') {
        central_tendency_values <- apply(M_ij, 1, mean)
      } else if (central_tendency == 'median') {
        central_tendency_values <- apply(M_ij, 1, median)
      } else {
        stop('Invalid central_tendency argument. Choose "mean" or "median".')
      }

      # Compute the lower and upper quantiles for the confidence bands
      lower_quantile <- significance_level / 2
      upper_quantile <- 1 - (significance_level / 2)

      lower_bounds <- apply(M_ij, 1, quantile, probs = lower_quantile)
      upper_bounds <- apply(M_ij, 1, quantile, probs = upper_quantile)

      # Create a data frame for plotting
      df <- data.frame(
        time = 1:T,
        central_tendency = central_tendency_values,
        lower_bound = lower_bounds,
        upper_bound = upper_bounds
      )

      # Create the plot
      p <- ggplot(df, aes(x = time)) +
        geom_line(aes(y = central_tendency)) +
        geom_line(aes(y = lower_bound), color = 'blue', linetype = 'dotted') +
        geom_line(aes(y = upper_bound), color = 'blue', linetype = 'dotted') +
        geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = 'pink', alpha = 0.3) +
        xlab('')+ylab('')
      # if(i==1){
        p<-p+labs(title =paste0('Shock: ',shock_names[j]),subtitle = paste0('Variable: ',variable_names[i]))
      # }
      # if(j==1){
        # p<-p+ylab(variable_names[i])
      # }
      p<-p +
        theme(title = element_text(size = 8, family = 'mono')) +
        geom_hline(yintercept = 0, color = 'red', linetype = 'dashed')+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))

      # Add the plot to the list
      plot_list[[plot_counter]] <- p
      plot_counter <- plot_counter + 1
    }
  }

  # Arrange the plots in a grid
  # Calculate the number of rows and columns based on the number of variables and shocks
    # grid_arranged_plots <- arrangeGrob(grobs = plot_list, nrow = N, ncol = s)

  # Return the list of plots
  return(plot_list)
}
