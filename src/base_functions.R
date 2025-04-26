#' Create data matrices for prices and quantities
#'
#' @param data The input data frame.
#' @param period_var The period variable.
#' @param variety_var The variety variable.
#' @param price_var The price variable.
#' @param quantity_var The quantity variable.
#' @return A list containing the price and quantity matrices.
create_data_matrices <- function(data, period_var, variety_var, price_var, quantity_var) {
  # Create price and quantity matrices for each period and variety
  price_matrix <- tapply(data[[price_var]], list(data[[period_var]], data[[variety_var]]), identity)
  quantity_matrix <- tapply(data[[quantity_var]], list(data[[period_var]], data[[variety_var]]), identity)
  
  return(list(price_matrix = price_matrix,
              quantity_matrix = quantity_matrix))
}

#' Calculate Tornqvist matrix
#'
#' @param price_matrix The price matrix.
#' @param quantity_matrix The quantity matrix.
#' @param window The window size for the calculation.
#' @return A matrix containing the Tornqvist indices.
calculate_tornqvist_matrix <- function(price_matrix, quantity_matrix, window) {
  # Calculate sales matrix
  sales_matrix <- price_matrix * quantity_matrix
  length_periods <- nrow(sales_matrix)
  
  # Initialize Tornqvist matrix to hold bilateral indices
  tornqvist_matrix <- matrix(NA_real_, length_periods, length_periods)
  diag(tornqvist_matrix) <- 0
  
  # Calculate bilateral Tornqvist indices
  for (r in 1:(length_periods - 1)) {
    for (c in (r + 1):min(r + window - 1, length_periods)) {
      # Calculate Tornqvist index for matched varieties [could be Fisher index]
      tornqvist_index <- sum(0.5 * (
        (sales_matrix[r, ] / sum(sales_matrix[r, ] * !is.na(sales_matrix[c, ]), na.rm = TRUE)) +
        (sales_matrix[c, ] / sum(sales_matrix[c, ] * !is.na(sales_matrix[r, ]), na.rm = TRUE))
      ) * log(price_matrix[c, ] / price_matrix[r, ]), na.rm = TRUE)
      
      # Populate matrix symmetrically using time reversibility
      tornqvist_matrix[r, c] <- tornqvist_index
      tornqvist_matrix[c, r] <- -tornqvist_index
    }
  }
  
  return(tornqvist_matrix)
}

#' Calculate GEKS index
#'
#' @param tornqvist_matrix The Tornqvist matrix.
#' @param window The window size for the calculation.
#' @return A numeric vector containing the GEKS index.
calculate_geks_index <- function(tornqvist_matrix, window) {
  # Initialize GEKS index vector
  length_periods <- nrow(tornqvist_matrix)
  geks_index <- numeric(length_periods)
  
  # Calculate GEKS index and apply mean splice
  for (t in 1:(length_periods - window + 1)) {
    # Calculate GEKS index for the current vintage
    geks_vintage <- colMeans(tornqvist_matrix[t:(t + window - 1), t:(t + window - 1)])
    
    # Apply mean splice to combine GEKS indices
    if (t == 1) {
      # For the first vintage, set GEKS index relative to the first period
      geks_index[t:(t + window - 1)] <- geks_vintage - geks_vintage[1]
    } else {
      # For subsequent vintages, extend GEKS index using mean splice [could be half splice]
      geks_index[t + window - 1] <- mean(geks_index[t:(t + window - 2)]) +
        (geks_vintage[window] - mean(geks_vintage[1:(window - 1)]))
    }
  }
  
  # Normalize GEKS index
  return(100 * exp(geks_index))
}