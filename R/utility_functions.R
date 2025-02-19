#' Remove Outliers from Data
#'
#' @description
#' This function filters out outliers from a data frame with two columns (`x` and `y`)
#' using the interquartile range (IQR) method. Any data points outside 1.5 * IQR
#' from the first and third quartiles are considered outliers.
#' Not exported for user access.
#'
#' @param data A data frame with two numeric columns `x` and `y`.
#'
#' @return A data frame with outliers removed.
#'
#' @importFrom stats quantile
#'
#' @examples
#' data <- data.frame(x = rnorm(100), y = rnorm(100))
#' PopComm:::data_clean <- remove_outlier(data)
#'
#' @keywords internal
#' @noRd
remove_outlier <- function(data) {
  # Calculate IQR for x and y
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


#' Project Points onto a Line
#'
#' @description
#' This function projects a point `(x, y)` onto a line defined by the equation
#' `y = slope * x + intercept`, and returns the coordinates of the projection.
#' Not exported for user access.
#'
#' @param x A numeric value representing the x-coordinate of the point to project.
#' @param y A numeric value representing the y-coordinate of the point to project.
#' @param slope A numeric value representing the slope of the line.
#' @param intercept A numeric value representing the intercept of the line.
#'
#' @return A numeric vector with the projected `(x, y)` coordinates on the line.
#'
#' @examples
#' PopComm:::project_to_line(1, 2, 3, 4)
#'
#' @keywords internal
#' @noRd
project_to_line <- function(x, y, slope, intercept) {
  # Calculate projection of point onto the line
  x_proj <- (x + slope * y - slope * intercept) / (slope^2 + 1)
  x_proj <- round(x_proj, 5)
  y_proj <- slope * x_proj + intercept
  y_proj <- round(y_proj, 5)

  return(c(x_proj, y_proj))
}
