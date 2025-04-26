#' Calculate M-GEKS index from prices and quantities
#'
#' @param data The input data frame.
#' @param period_var The period variable.
#' @param variety_var The variety variable.
#' @param price_var The price variable.
#' @param quantity_var The quantity variable.
#' @param window The window size for the calculation.
#' @param use_rcpp A logical value indicating whether to use Rcpp.
#' @return A data frame containing the GEKS index for each period.
matrix_geks <- function(data, period_var, variety_var, price_var, quantity_var, window, use_rcpp = FALSE) {
  # Convert variable names to strings
  period_var <- deparse(substitute(period_var))
  variety_var <- deparse(substitute(variety_var))
  price_var <- deparse(substitute(price_var))
  quantity_var <- deparse(substitute(quantity_var))
  
  # Sort unique periods
  periods <- sort(unique(data[[period_var]]))
  
  # Call functions using Rcpp or base R
  if (use_rcpp) {
    data_matrices <- create_data_matrices_cpp(data, period_var, variety_var, price_var, quantity_var)
    tornqvist_matrix <- calculate_tornqvist_matrix_cpp(data_matrices$price_matrix, data_matrices$quantity_matrix, window)
    geks_index <- calculate_geks_index_cpp(tornqvist_matrix, window)
  } else {
    data_matrices <- create_data_matrices(data, period_var, variety_var, price_var, quantity_var)
    tornqvist_matrix <- calculate_tornqvist_matrix(data_matrices$price_matrix, data_matrices$quantity_matrix, window)
    geks_index <- calculate_geks_index(tornqvist_matrix, window)
  }
  
  return(data.frame(period = periods,
                    index = geks_index))
}