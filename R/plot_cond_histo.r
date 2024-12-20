#' plot_cond_histo function
#'
#' This functiion uses the conditional probability calculations (eg scenarios) and plots the histogram of the selected variable
#'
#' @param variable Name of variable to be plotted (string)
#' @param horizon At which horizon (horizon<=h)
#' @param threshold (optional) If present compute P(x>threshold)
#'
#' @returns ggplot object (plot)
#'
#' @export
#'
#' @import dplyr

plot_cond_histo<-function(variable=NULL,horizon=1,threshold=NULL){
  y_h_df<-as.data.frame(t(y_h))
  comby=paste(variable,horizon,sep='.')

  p <- ggplot(data = y_h_df) +
    geom_histogram(aes(x = !!sym(comby)), alpha = 0.5) +
    labs(title = paste0('Distribution of forecast at horizon ', horizon)) +
    theme_minimal()

  if (!is.null(threshold)) {

    # Calculate the probability of values above the threshold
    med_=median(y_h_df[[comby]])
    if(threshold>med_){
    prob_above_threshold <- mean(y_h_df[[comby]] > threshold)
    }else{
      prob_above_threshold <- mean(y_h_df[[comby]] < threshold)
    }
    # Add the vertical line and annotate the probability on the plot
    p <- p +
      geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
      annotate("text", x = threshold, y = 0, label = paste0("P(X > ", threshold, ") = ", round(prob_above_threshold, 3)),
               vjust = -0.5, hjust = -0.1, color = "red")
  }

  return(p)

  }
