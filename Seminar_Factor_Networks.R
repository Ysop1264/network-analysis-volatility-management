##### Uncomment if packages not installed
install.packages("tidyverse")
install.packages("tidyfinance")
install.packages("scales")
install.packages("frenchdata")
install.packages("dplyr")
install.packages("moments")
install.packages("sandwich")
install.packages("lmtest")
install.packages("lubridate")
install.packages("nlshrink")
install.packages("rugarch")
install.packages("xts")
install.packages("zoo")
#####
library(rlang)
library(tidyverse)
library(tidyfinance)
library(scales)
library(frenchdata)
library(dplyr)
library(moments)
library(sandwich)
library(lmtest)
library(lubridate)
library(nlshrink)
library(rugarch)
library(xts)
library(zoo)

# DATA download and transformation
start_date <- as.Date("1971-01-01")
end_date   <- as.Date("2026-01-31")

# Code from tidy-fianance website to download data
factors_ff5_daily <- download_french_data("Fama/French 5 Factors (2x3) [Daily]")$subsets$data[[1]] |>
  mutate(
    date = ymd(date), 
    across(c(RF, "Mkt-RF", SMB, HML, RMW, CMA), ~ as.numeric(.)/100), 
    .keep = "none"
  ) |> 
  rename_with(str_to_lower) |>
  rename(mkt_excess = "mkt-rf", risk_free = rf) |>
  filter(date >= start_date & date <= end_date)

momentum_daily <- download_french_data("Momentum Factor (Mom) [Daily]")$subsets$data[[1]] |>
  mutate(
    date = ymd(date),
    across(Mom, ~ as.numeric(.)/100),
    .keep = "none"
  ) |>
  rename_with(str_to_lower) |>
  filter(date >= start_date & date <= end_date)

industry_daily <- download_french_data("49 Industry Portfolios [Daily]")$subsets$data[[1]] |>
  mutate(
    date = ymd(date), 
    across(-date, ~ as.numeric(.)/100), 
    .keep = "unused"
  ) |> 
  rename_with(str_to_lower) |>
  filter(date >= start_date & date <= end_date)

factors_joined <- factors_ff5_daily |> left_join(momentum_daily, by = "date") |>
  left_join(industry_daily, by = "date")

factors_joined_excess <- factors_joined |> 
  mutate(
    across(-c(mkt_excess, smb, hml, rmw, cma, risk_free), ~ . - risk_free)
  )

managed_portfolios <- factors_joined_excess |> select(-risk_free)

print(head(managed_portfolios))

# =================================================
# Functions
# =================================================

#' Compute the MVE static weights
#' 
#' 
#' @param mu Vector of means
#' @param sigma Variance-covariance matrix of returns
#' @return Solved weights b
compute_MVE_weights <- function(mu, sigma){
  b <- solve(sigma, mu)/ sqrt(as.numeric(t(mu)%*%solve(sigma)%*%mu))
  return (b)
}

#' Computes the Sharpe Ratio of the returns
#' 
#' 
#' @param r Vector of returns
#' @return Sharpe Ratio
compute_SR <- function(r){
  SR = mean(r)/sd(r)
  return(SR)
}

# Constructed Realised Variance
rv_monthly <- managed_portfolios %>%
  mutate(month = floor_date(date, "month")) %>%
  pivot_longer(
    cols = -c(date, month),
    names_to = "factor",
    values_to = "ret"
  ) %>%
  group_by(factor, month) %>%
  summarise(
    RV = sum((ret - mean(ret, na.rm = TRUE))^2, na.rm = TRUE),
    sigma = sqrt(RV),
    .groups = "drop"
  )

print(head(rv_monthly))
print(tail(rv_monthly))

# DCC-NL Estimation for Covariance Matrix
# Fitting univariate GARCH(1,1) 
# Function to fit GARCH(1,1) for one return series on an expanding window (daily estimation)
fit_garch_expanding <- function(
  ret, dates = NULL, 
  init_window = 504, # 2 trading years approx
  refit_on = c("month_end", "every_day"),
  garch_order = c(1,1),
  arma_order = c(0,0),
  distribution = "norm"
  ) {
  # Handling inputs
  ret <- as.numeric(ret)
  dates <- as.Date(dates)

  n <- length(ret)
  if (n <= init_window) {
    stop("Series shorter than initial window.")
  }

  # Specifying GARCH model
  spec <- ugarchspec(
    variance.model = list(
      model = "sGARCH", garchOrder = garch_order
    ), mean.model = list(
      armaOrder = arma_order,
      include.mean = TRUE
    ), 
    distribution.model = distribution
  )

  # Storing results
  results <- vector("list", n - init_window)

  # Expanding estimation
  for (t in (init_window):(n - 1)){
    # The set of insample returns
    insample_ret <- ret[1:t]

    # Fitting GARCH based on insample data
    fit <- ugarchfit(spec = spec, data = insample_ret, solver = "Hybrid")

    # To catch errora
    if (is.null(fit)) {
      results[[t - init_window + 1]] <- data.frame(
        date = dates[t + 1],
        mu = NA_real_,
        sigma = NA_real_,
        sigma2 = NA_real_,
        resid = NA_real_,
        z = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
  }

  # 1-step ahead forecast
  fc <- ugarchforecast(fit, n.ahead = 1)
  mu_fc <- as.numeric(fitted(fc))
  sigma_fc <- as.numeric(sigma(fc))
  sigma2_fc <- sigma_fc^2

  # Realised next return
  r_next <- ret[t + 1]

  # standardized residual (for DCC later)
    z_next <- (r_next - mu_fc) / sigma_fc

    results[[t - init_window + 1]] <- data.frame(
      date = dates[t + 1],
      mu = mu_fc,
      sigma = sigma_fc,
      sigma2 = sigma2_fc,
      resid = r_next - mu_fc,
      z = z_next,
      stringsAsFactors = FALSE
    )

  out <- do.call(rbind, results)

  return(out)
}

# Adjacency matrix construction 
# @param sigma_hat is the forecasted covariance matrix 
# @param tau is the threshold (0.1,0.5)
# @return List with partial correlations and the adjacency matrix
adjacency_matrix <- function(sigma_hat, tau = 0.1) {
  # start with precision matrix 
  theta <- solve(sigma_hat)
  N <- nrow(theta)
  # vector approach to get the partial correlation matrix
  # make a vector for the inverse of the square roots of the diagonal elements
  d <- 1/sqrt(diag(theta))
  
  # partial correlation is -theta_ij / sqrt(theta_ii * theta_jj)
  partial_corr <- -theta * (d %o% d)
  # set diagonals to 0 to prevent them from being -1 with themselves
  diag(partial_corr) <- 0

  # now separate positive and negative adjacency matrix
  # positive (contagion):
  adjacency_pos <- ifelse(partial_corr > tau, partial_corr, 0)
  # negative (potential for hedging):
  adjacency_neg <- ifelse(partial_corr < -tau, -partial_corr, 0)
  return(list(
    partial_corr = partial_corr,
    adjacency_pos = adjacency_pos,
    adjacency_neg = adjacency_neg
  ))
}

# Grid search for tau
# define the grid 
taus <- seq(0.1, 0.5, by = 0.1)
# store the results for each tau 
grid_tau <- list()
for (i in taus) {
  # execute the function for the current tau starting from 0.1
  adjacency_output <- adjacency_matrix(sigma_hat, tau = i)
  grid_tau[[as.character(i)]] <- adjacency_output
  count_positive <- sum(adjacency_output$adjacency_pos > 0)
  count_negative <- sum(adjacency_output$adjacency_neg > 0)

  # check
  cat(sprintf("Tau: %.1f | Pos Links: %d | Neg Links: %d\n", i, count_positive, count_negative))
}

# Network Centrality function using EC 
# @param adjacency_pos for positive adjacency matrix
# @param adjacency_neg for negative adjacency matrix 
# @return List with EC+ and EC-
network_centrality <- function(adjacency_pos, adjacency_neg){
  # eigenvalue for positive (contagion) network:
  eigenvalue_pos <- eigen(adjacency_pos)
  # eigenvector corresponding to largest eigenvalue for positive (contagion):
  eigenvector_pos <- Re(eigenvalue_pos$vectors[,1])
  # ensure positive values 
  if (sum(eigenvector_pos) < 0) eigenvector_pos <- -eigenvector_pos
  # normalize them 
  ec_pos <- eigenvector_pos / sum(eigenvector_pos)
  
  # eigenvalue for negative (potential for hedging) network:
  eigenvalue_neg <- eigen(adjacency_neg)
  # eigenvector corresponding to largest eigenvalue for negative (hedging):
  eigenvector_neg <- Re(eigenvalue_neg$vectors[,1])
  # ensure positive values 
  if (sum(eigenvector_neg) < 0) eigenvector_neg <- -eigenvector_neg
  # normalize them
  ec_neg <- eigenvector_neg / sum(eigenvector_neg)

  #return list 
  return(list(
    ec_pos = ec_pos,
    ec_neg = ec_neg
  ))
}



# =================================================
# Benchmarks
# =================================================

# Setting here to test the script for less managed portfolios
managed_portfolios <- managed_portfolios[,1:7]

# Defining the estimation sample
start_date_estimation <- as.Date("1971-01-01")
end_date_estimation <- as.Date("1973-01-01")


### MVE Benchmark by Moreira and Muir (2017)

# Calculating in-sample mean
managed_portfolios_est <- managed_portfolios |> 
  filter(date >= start_date_estimation & date <= end_date_estimation)
mu <- as.matrix(colMeans(managed_portfolios_est[,-1]))
colnames(mu) <- c("mean")

# Calculating the in-sample covariance matrix with Shrinkage estimator (Ledoit and Wolf (2004))
sigma <- as.matrix(linshrink_cov(as.matrix(managed_portfolios_est[,-1])))

# Calculating the optimal weights 
b <- compute_MVE_weights(mu, sigma)

# Resulting in-sample maximized sharpe ratio 
sp_in_sample <- t(b)%*%mu
print(paste("In-sample Sharpe Ratio:",  round(sp_in_sample,4)))


# The resulting weighted factor before volatility management 
MVE_returns <- as.matrix(managed_portfolios[,-1]) %*% b
MVE_returns_df <- tibble(
  date = managed_portfolios$date,
  month = floor_date(date, "month"),
  MVE_return = as.numeric(MVE_returns)
  
)

# Computing the realized variance
rv_lag_MVE <- MVE_returns_df |>
  group_by(month) |>
  summarise(
    n_of_days = n(),
    rv = sum((MVE_return-mean(MVE_return))^2) / n_of_days,
    .groups = "drop"
  ) |> mutate(
    rv_lag = lag(rv)
  ) |> select(-rv)

# Scaling the portfolio by its lagged realized variance (REPAIR THAT BC I NEED TO SCALE BY LAG!!!)
MVE_returns_df <- MVE_returns_df |> left_join(rv_lag_MVE, by = "month")
MVE_returns_df <- MVE_returns_df |> 
  mutate(
    MVE_scaled = MVE_return / rv_lag
  )
MVE_returns_df <- MVE_returns_df |> select(-c("rv_lag", "n_of_days"))

# Take out the first month, because we do not have the realized variance for the period 0
MVE_returns_df <- MVE_returns_df |> filter(month >= start_date_estimation + months(1))

# Scaling by c
c = sd(MVE_returns_df$MVE_return)/sd(MVE_returns_df$MVE_scaled)
MVE_returns_df$MVE_strategy_return <- MVE_returns_df$MVE_scaled*c
MVE_returns_df <- MVE_returns_df |> select(-c(MVE_return, MVE_scaled))


### Equally weighted buy-and-hold 
no_of_factors <- dim(managed_portfolios[,-1])[2]
b_EW <- as.matrix(rep(1/no_of_factors, no_of_factors))

EW_returns <- as.matrix(managed_portfolios[,-1]) %*% b_EW
EW_returns_df <- tibble(
  date = managed_portfolios$date, 
  month = floor_date(date, "month"),
  EW_return = as.numeric(EW_returns)
) |> filter(month >= start_date_estimation + months(1))


# Changing data frames only for out of sample periods
MVE_returns_df <- MVE_returns_df |> filter(month > end_date_estimation)
EW_returns_df <- EW_returns_df |> filter(month > end_date_estimation)



## This is just for veryfying the performance of the benchmarks

# Calculating the Sharpe Ratios
sharpe_ratios_benchmarks <- data_frame(
  Benchmark = c("EW Buy-and-Hold", "MVE"),
  SR = c(round(compute_SR(EW_returns_df$EW_return),4)*sqrt(252), 
         round(compute_SR(MVE_returns_df$MVE_strategy_return),4)*sqrt(252))
)


# Creating data frame for letter regressions
benchmarks_returns <- EW_returns_df |> left_join(MVE_returns_df, by = "date") |> select(-c(month.x, month.y))
