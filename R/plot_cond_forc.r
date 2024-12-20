#' plot_cond_forc function
#'
#' @param variable name of variable to be plotted (string)
#'
#' @returns plot
#'
#' @export
#'
#' @import dplyr
#'
plot_cond_forc<-function(variable=varbls[1]){
  p <- ggplot(cond.for[cond.for$variable == variable, ], aes(x = hor)) +
    # Median line (solid line)
    geom_line(aes(y = center, color =factor(hist),group=hist), linewidth = 1, show.legend = TRUE) +
    # Shaded area for 68% HDI
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "bnd"),
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

    # Custom legend for fill
    #scale_fill_manual(
    #  name = "Legend",
    #  labels=c('bnd'="68% HDI (Highest Density Interval) total"),
    #  values = c('bnd' = "lightblue")
    #)+
    theme(legend.position = 'none')


  # Display the plot

  return(p)
}
