##### Uncomment if packages not installed
# install.packages("tidyverse")
# install.apckages("tidyfinance")
# install.packages("scales")
# install.packages("frenchdata")
# install.packages("dplyr")
# install.packages("moments")
# install.packages("sandwich")
# install.packages("lmtest")
# install.packages("lubridate")
# install.packages("nlshrink")

#####

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
