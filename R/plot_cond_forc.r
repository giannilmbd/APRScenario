#' plot_cond_forc function; Data should conatain the variable "variable", the "hor" horizon and a "history"
#'
#' @param variable name of variable to be plotted (string)
#' @param data: conditional forecast data frame (eg cond.for)
#' @returns plot
#'
#' @export
#'
#' @import dplyr
#'
plot_cond_forc<-function(variable=NULL,data=NULL){
  cond.for<-data
  p <- ggplot(cond.for[cond.for$variable == variable, ], aes(x = hor)) +
    # Median line (solid line)
    geom_line(aes(y = center, color =factor(hist),group=hist), linewidth = 1, show.legend = TRUE) +
    # Shaded area for 68% HDI
    geom_ribbon(aes(ymin = lower, ymax = upper),fill='pink',
                alpha = 0.5, show.legend = TRUE) +
    # Labels and theme
    labs(title = "Conditional (scenario) Forecast", x = "h", y = variable) +
    theme_minimal() +
    # Custom legend for colors
    scale_color_manual(
      name = "",
      labels=c('0'="",'1'=""),
      values = c("0" = "blue",  '1'= "red")
    ) +
    theme(legend.position = 'none')


  # Display the plot

  return(p)
}
