# utilities
remove_outlier <- function(data) {
  Q1_x <- quantile(data$x, 0.25, na.rm = TRUE)
  Q3_x <- quantile(data$x, 0.75, na.rm = TRUE)
  IQR_x <- Q3_x - Q1_x
  
  Q1_y <- quantile(data$y, 0.25, na.rm = TRUE)
  Q3_y <- quantile(data$y, 0.75, na.rm = TRUE)
  IQR_y <- Q3_y - Q1_y
  
  # Define outlier thresholds
  lower_bound_x <- Q1_x - 1.5 * IQR_x
  upper_bound_x <- Q3_x + 1.5 * IQR_x
  lower_bound_y <- Q1_y - 1.5 * IQR_y
  upper_bound_y <- Q3_y + 1.5 * IQR_y
  
  # Identify non-outliers
  non_outliers <- (data$x > lower_bound_x) & (data$x < upper_bound_x) &
    (data$y > lower_bound_y) & (data$y < upper_bound_y)

  # Filter data to remove outliers
  data_clean <- data[non_outliers, ]
  return(data_clean)
}


project_to_line <- function(x, y, slope, intercept) {
  x_proj <- (x + slope * y - slope * intercept) / (slope^2 + 1)
  y_proj <- slope * x_proj + intercept
  return(c(x_proj, y_proj))
}
