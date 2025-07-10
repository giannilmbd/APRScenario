#' plot_cond_histo function
#'
#' This function uses the conditional probability calculations (eg scenarios) and plots the histogram of the selected variable
#'
#' @param data data of conditional forecasts
#' @param variable (character) Name of variable to be plotted
#' @param horizon (numeric) At which horizon (horizon<=h)
#' @param threshold (numeric,optional) If present compute P(x>threshold)
#' @param above (logical,optional): if TRUE then compute probability above threshold
#' @param size (optional) size of annotation text in the plot
#'
#' @returns ggplot object (plot)
#' @examples
#' # Example with simulated conditional forecast data
#' # Create sample forecast data matrix
#' set.seed(123)
#' n_sims <- 500
#' horizons <- 3
#' variables <- c("GDP", "CPI", "FFR")
#' 
#' # Create column names in the expected format (variable.horizon)
#' col_names <- outer(variables, 1:horizons, paste, sep = ".")
#' 
#' # Generate random forecast data
#' forecast_data <- matrix(rnorm(n_sims * length(col_names)), 
#'                        nrow = n_sims, ncol = length(col_names))
#' colnames(forecast_data) <- as.vector(col_names)
#' 
#' # Plot histogram for GDP at horizon 2
#' p <- plot_cond_histo(data = t(forecast_data), 
#'                      variable = "GDP", 
#'                      horizon = 2,
#'                      threshold = 0.5, 
#'                      above = TRUE)
#'
#' @export
#'
#' @import dplyr
#' @importFrom ggplot2 ggplot aes geom_histogram after_stat geom_density labs theme_minimal geom_vline annotate

plot_cond_histo<-function(variable=NULL,horizon=1,threshold=NULL,data=NULL,above=TRUE,size=5){
  y_h=data
  y_h_df<-as.data.frame(t(y_h))
  comby=paste(variable,horizon,sep='.')

  p <- ggplot(data = y_h_df,aes(x = !!sym(comby))) +
    geom_histogram(aes(y=after_stat(density)), alpha = 0.1,fill='cyan') +
    geom_density(color='blue',fill=NA)+
    labs(title = paste0('Distribution of forecast at horizon ', horizon)) +
    theme_minimal()
  if (is.null(above)) {above=TRUE}
  if (!is.null(threshold)) {

    # Calculate the probability of values above the threshold
    # med_=median(y_h_df[[comby]])
    if(above){
    prob_above_threshold <- mean(y_h_df[[comby]] > threshold)
    }else{
      prob_above_threshold <- mean(y_h_df[[comby]] < threshold)
    }
    # Add the vertical line and annotate the probability on the plot
    p <- p +
      geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
      if(above){
      annotate("text", x = threshold, y = 0, label = paste0("P(X > ", round(threshold,3), ") = ", round(prob_above_threshold, 3)),
               vjust = -0.5, hjust = -0.1, color = "red",size=size)
      }else{
        annotate("text", x = threshold, y = 0, label = paste0("P(X < ", round(threshold,3), ") = ", round(prob_above_threshold, 3)),
                 vjust = -0.5, hjust = -0.1, color = "red",size=size)
      }
  }

  return(p)

  }
