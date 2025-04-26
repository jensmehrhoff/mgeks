#include <RcppArmadillo/Lighter>
// [[Rcpp::depends(RcppArmadillo)]]

//' Create data matrices for prices and quantities
//'
//' @param data The input data frame.
//' @param period_var The period variable.
//' @param variety_var The variety variable.
//' @param price_var The price variable.
//' @param quantity_var The quantity variable.
//' @return A list containing the price and quantity matrices.
// [[Rcpp::export]]
Rcpp::List create_data_matrices_cpp(Rcpp::DataFrame data, std::string period_var, std::string variety_var, std::string price_var, std::string quantity_var) {
  // Extract columns from the data frame
  Rcpp::CharacterVector period = data[period_var];
  Rcpp::CharacterVector variety = data[variety_var];
  arma::colvec price = data[price_var];
  arma::colvec quantity = data[quantity_var];
  
  // Get unique periods and varieties
  Rcpp::CharacterVector unique_periods = sort_unique(period);
  Rcpp::CharacterVector unique_varieties = unique(variety);
  
  // Create maps to store period and variety indices
  std::unordered_map<std::string, int> period_map;
  std::unordered_map<std::string, int> variety_map;
  for (int i = 0; i < unique_periods.size(); i++) {
    period_map[std::string(unique_periods[i])] = i;
  }
  for (int i = 0; i < unique_varieties.size(); i++) {
    variety_map[std::string(unique_varieties[i])] = i;
  }
  
  // Initialize price and quantity matrices
  arma::mat price_matrix(unique_periods.size(), unique_varieties.size(), arma::fill::value(NA_REAL));
  arma::mat quantity_matrix(unique_periods.size(), unique_varieties.size(), arma::fill::value(NA_REAL));
  
  // Populate the matrices for each period and variety
  for (int i = 0; i < data.nrows(); i++) {
    int period_index = period_map[std::string(period[i])];
    int variety_index = variety_map[std::string(variety[i])];
    price_matrix(period_index, variety_index) = price[i];
    quantity_matrix(period_index, variety_index) = quantity[i];
  }
  
  return Rcpp::List::create(Rcpp::Named("price_matrix") = price_matrix,
                            Rcpp::Named("quantity_matrix") = quantity_matrix);
}

//' Calculate Tornqvist matrix
//'
//' @param price_matrix The price matrix.
//' @param quantity_matrix The quantity matrix.
//' @param window The window size for the calculation.
//' @return A matrix containing the Tornqvist indices.
// [[Rcpp::export]]
arma::mat calculate_tornqvist_matrix_cpp(const arma::mat& price_matrix, const arma::mat& quantity_matrix, int window) {
  // Calculate sales matrix
  arma::mat sales_matrix = price_matrix % quantity_matrix;
  int n_periods = sales_matrix.n_rows;
  int n_varieties = sales_matrix.n_cols;
  
  // Initialize Tornqvist matrix to hold bilateral indices
  arma::mat tornqvist_matrix(n_periods, n_periods, arma::fill::value(NA_REAL));
  tornqvist_matrix.diag().zeros();
  
  // Calculate bilateral Tornqvist indices
  for(int r = 0; r < n_periods - 1; ++r) {
    for(int c = r + 1; c < std::min(r + window, n_periods); ++c) {
      // Sum sales for periods r(ow) and c(olumn)
      double sum_sales_r = 0.0;
      double sum_sales_c = 0.0;
      for(int i = 0; i < n_varieties; ++i) {
        if(!std::isnan(sales_matrix(r, i)) && !std::isnan(sales_matrix(c, i))) {
          sum_sales_r += sales_matrix(r, i);
          sum_sales_c += sales_matrix(c, i);
        }
      }
      
      // Calculate Tornqvist index for matched varieties [could be Fisher index]
      double tornqvist_index = 0.0;
      for(int i = 0; i < n_varieties; ++i) {
        if(!std::isnan(sales_matrix(r, i)) && !std::isnan(sales_matrix(c, i))) {
          tornqvist_index += 0.5 * (
            (sales_matrix(r, i) / sum_sales_r) +
            (sales_matrix(c, i) / sum_sales_c)
          ) * std::log(price_matrix(c, i) / price_matrix(r, i));
        }
      }
      
      // Populate matrix symmetrically using time reversibility
      tornqvist_matrix(r, c) = tornqvist_index;
      tornqvist_matrix(c, r) = -tornqvist_index;
    }
  }
  
  return tornqvist_matrix;
}

//' Calculate GEKS index
//
//' @param tornqvist_matrix The Tornqvist matrix.
//' @param window The window size for the calculation.
//' @return A numeric vector containing the GEKS index.
// [[Rcpp::export]]
arma::colvec calculate_geks_index_cpp(const arma::mat& tornqvist_matrix, int window) {
  // Initialize GEKS index vector
  int n_periods = tornqvist_matrix.n_rows;
  arma::colvec geks_index(n_periods, arma::fill::value(NA_REAL));
  
  // Calculate GEKS index and apply mean splice
  for(int t = 0; t < n_periods - window + 1; ++t) {
    // Calculate GEKS index for the current vintage
    arma::colvec geks_vintage = arma::mean(tornqvist_matrix.submat(t, t, t + window - 1, t + window - 1), 0).t();
    
    // Apply mean splice to combine GEKS indices
    if (t == 0) {
      // For the first vintage, set GEKS index relative to the first period
      geks_index.subvec(0, window - 1) = geks_vintage - geks_vintage(0);
    } else {
      // For subsequent vintages, extend GEKS index using mean splice [could be half splice]
      geks_index(t + window - 1) = arma::mean(geks_index.subvec(t, t + window - 2)) +
        (geks_vintage(window - 1) - arma::mean(geks_vintage.subvec(0, window - 2)));
    }
  }
  
  // Normalize GEKS index
  return 100 * arma::exp(geks_index);
}