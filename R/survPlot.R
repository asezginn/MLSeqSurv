survPlot <- function(model){

  survHazards <- cbind(rep(1, nrow(model@survHazards)), model@survHazards)
  event_times <- model@times
  event_times <- c(0, event_times)

  plot_list <- list()

  for (i in 1:nrow(survHazards)) {
    survival_curve <- survHazards[i, ]

    p <- ggplot(data = data.frame(event_times = event_times, survival_curve = survival_curve),
                aes(x = event_times, y = survival_curve)) +
      geom_step() +
      labs(x = "Time", y = "Survival probability", title = paste("Survival Curve for Person", i)) +
      ylim(0, 1)

    plot_list[[i]] <- p
  }


  grid.arrange(grobs = plot_list, ncol = 4)

}
