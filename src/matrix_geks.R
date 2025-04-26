#' A Rapid and Efficient Compilation Procedure for the GEKS Index: The Matrix-Based (M-GEKS) Approach
#' https://unece.org/info/events/event/394156
#'
#' (C) 2025 by Jens Mehrhoff
#' jens.mehrhoff(at)bundesbank.de
#'
#' Licensed under CC BY 4.0
#' https://creativecommons.org/licenses/by/4.0/
#' 
#' Analysis-ready data exported from Dominick's Finer Foods data set, https://github.com/eurostat/dff
#' Jens Mehrhoff (2018), Promoting the Use of a Publicly Available Scanner Data Set

# SET-UP R ----------------------------------------------------------------

# Install packages and set options
load_packages <- function(packages) {
  for (package in packages) {
    if(!require(package, character.only = TRUE)) {
      install.packages(package)
      library(package, character.only = TRUE)
    }
  }
}
packages <- c("tidyverse", "Rcpp", "RcppArmadillo", "IndexNumR", "PriceIndices", "bench")
load_packages(packages)
options(warn = 0)
rm(list = ls())

# Set working directory
setwd("/cloud/project")

# LOAD DATA ---------------------------------------------------------------

# Read analysis-ready data
move_monthly <- read_csv("ard.csv") %>%
  filter(between(MONTH, ymd("1990-01-01"), ymd("1997-04-01")))

# Calculate unit prices and create product codes
move_monthly <- move_monthly %>%
  group_by(CAT) %>%
  mutate(
    PRICE = SALES / MOVE,
    PRODUCT = interaction(COM_CODE, NITEM, drop = TRUE) %>% as.character()
  ) %>%
  ungroup()

# Introduce artificial final period
move_monthly <- move_monthly %>%
  bind_rows(
    move_monthly %>%
      filter(MONTH == ymd("1990-01-01")) %>%
      mutate(MONTH = ymd("1997-05-01"))
  )

# CALCULATE PRICE INDICES ---------------------------------------------------------

# Source main function as well as base R and Rcpp functions
source("main_function.R")
source("base_functions.R")
Rcpp::sourceCpp("rcpp_functions.cpp")

# Base R
system.time(
  geks_index_base <- move_monthly %>%
    split(.$CAT) %>%
    map(~ matrix_geks(
      data = .x,
      period_var = MONTH,
      variety_var = PRODUCT,
      price_var = PRICE,
      quantity_var = MOVE,
      window = 25L,
      use_rcpp = FALSE
    )) %>%
    list_rbind(names_to = "category")
)

# Rcpp
system.time(
  geks_index_rcpp <- move_monthly %>%
    split(.$CAT) %>%
    map(~ matrix_geks(
      data = .x,
      period_var = MONTH,
      variety_var = PRODUCT,
      price_var = PRICE,
      quantity_var = MOVE,
      window = 25L,
      use_rcpp = TRUE
    )) %>%
    list_rbind(names_to = "category")
)

stopifnot(all.equal(geks_index_base, geks_index_rcpp))

# IndexNumR
system.time(
  geks_index_indexnumr <- move_monthly %>%
    mutate(
      across(MONTH, ~ interval(ymd("1990-01-01"), .x) %/% months(1) + 1)
    ) %>%
    split(.$CAT) %>%
    map(~ IndexNumR::GEKSIndex(
      x = .x,
      pvar = "PRICE",
      qvar = "MOVE",
      pervar = "MONTH",
      prodID = "PRODUCT",
      indexMethod = "tornqvist",
      window = 25L,
      splice = "mean_pub"
    ) %>%
      as.data.frame() %>%
      mutate(
        period = ymd("1990-01-01") + months(row_number() - 1),
        index = 100 * V1
      ) %>%
      select(period, index)) %>%
    list_rbind(names_to = "category")
)

stopifnot(all.equal(geks_index_base, geks_index_indexnumr))

# PriceIndices
system.time(
  geks_index_priceindices <- move_monthly %>%
    rename(
      time = MONTH,
      prodID = PRODUCT,
      prices = PRICE,
      quantities = MOVE
    ) %>%
    split(.$CAT) %>%
    map(~ PriceIndices::price_indices(
      data = .x,
      start = "1990-01",
      end = "1997-05",
      formula = "ccdi_splice",
      window = 25L,
      splice = "mean_published",
      interval = TRUE
    ) %>%
      mutate(
        period = ym(time),
        index = 100 * ccdi_splice
      ) %>%
      select(period, index)) %>%
    list_rbind(names_to = "category")
)

stopifnot(all.equal(geks_index_base, geks_index_priceindices))

# BENCHMARK INDEX FUNCTIONS -----------------------------------------------

# Base R (by category)
bench_mark_base <-move_monthly %>%
  split(.$CAT) %>%
  map(~ bench::mark(
    matrix_geks(
      data = .x,
      period_var = MONTH,
      variety_var = PRODUCT,
      price_var = PRICE,
      quantity_var = MOVE,
      window = 25L,
      use_rcpp = FALSE
    ),
    iterations = 100,
    check = FALSE,
    memory = FALSE,
    filter_gc = FALSE
  )) %>%
  list_rbind(names_to = "category")

# Rcpp (by category)
bench_mark_rcpp <-move_monthly %>%
  split(.$CAT) %>%
  map(~ bench::mark(
    matrix_geks(
      data = .x,
      period_var = MONTH,
      variety_var = PRODUCT,
      price_var = PRICE,
      quantity_var = MOVE,
      window = 25L,
      use_rcpp = TRUE
    ),
    iterations = 100,
    check = FALSE,
    memory = FALSE,
    filter_gc = FALSE
  )) %>%
  list_rbind(names_to = "category")
