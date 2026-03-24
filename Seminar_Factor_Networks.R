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
# install.packages("rugarch")
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
library(rugarch)

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
# Function to fit GARCH(1,1) for one return series on an expanding window (monthly estimation)
fit_garch_expanding_monthly <- function(
  ret,
  dates = NULL,
  init_window = 504,          # ~2 trading years
  garch_order = c(1, 1),
  arma_order = c(0, 0),
  distribution = "norm",
  scale_ret = TRUE
) {
  # Basic input handling
  ret <- as.numeric(ret)
  dates <- as.Date(dates)

  if (length(ret) != length(dates)) {
    stop("ret and dates must have the same length.")
  }

  if (any(is.na(dates))) {
    stop("dates contains NA values.")
  }

  if (scale_ret) {
    ret <- 100 * ret
  }

  n <- length(ret)

  if (n <= init_window) {
    stop("Series shorter than initial window.")
  }

  # GARCH specification
  spec <- ugarchspec(
    variance.model = list(
      model = "sGARCH",
      garchOrder = garch_order
    ),
    mean.model = list(
      armaOrder = arma_order,
      include.mean = TRUE
    ),
    distribution.model = distribution
  )

  # Safe wrappers
  # Added because of initial convergence issues
  safe_fit <- function(x) {
    tryCatch(
      suppressWarnings(
        ugarchfit(
          spec = spec,
          data = x,
          solver = "hybrid",
          solver.control = list(trace = 0)
        )
      ),
      error = function(e) NULL
    )
  }

  safe_forecast <- function(fit, n_ahead) {
    tryCatch(
      ugarchforecast(fit, n.ahead = n_ahead),
      error = function(e) NULL
    )
  }

  # Find monthly refit dates
  # Use last available trading day of each month
  # After init_window
  idx_all <- seq_len(n)
  month_id <- format(dates, "%Y-%m")

  # last observation index in each month
  month_end_idx <- tapply(idx_all, month_id, max)
  month_end_idx <- as.integer(month_end_idx)

  # only refit when we have enough history
  # and only if there is at least one next-day forecast to make
  refit_idx <- month_end_idx[month_end_idx >= init_window & month_end_idx < n]

  if (length(refit_idx) == 0) {
    stop("No eligible monthly refit dates found after init_window.")
  }

  # Prepare output container
  # Forecasts are for dates (init_window + 1):n
  out_list <- vector("list", n - init_window)
  out_pos <- 1

  # Loop over monthly refit dates
  for (k in seq_along(refit_idx)) {
    t_refit <- refit_idx[k]

    insample_ret <- ret[1:t_refit]

    # Determine how many future days this fit should cover
    next_refit <- if (k < length(refit_idx)) refit_idx[k + 1] else n

    # Forecast for days (t_refit + 1) up to next_refit
    forecast_idx <- (t_refit + 1):next_refit

    if (length(forecast_idx) == 0) next

    # Skip degenerate window
    if (sd(insample_ret, na.rm = TRUE) <= 1e-8) {
      for (j in forecast_idx) {
        out_list[[out_pos]] <- data.frame(
          date = dates[j],
          mu = NA_real_,
          sigma = NA_real_,
          sigma2 = NA_real_,
          resid = NA_real_,
          z = NA_real_,
          refit_date = dates[t_refit]
        )
        out_pos <- out_pos + 1
      }
      next
    }

    fit <- safe_fit(insample_ret)

    if (is.null(fit)) {
      for (j in forecast_idx) {
        out_list[[out_pos]] <- data.frame(
          date = dates[j],
          mu = NA_real_,
          sigma = NA_real_,
          sigma2 = NA_real_,
          resid = NA_real_,
          z = NA_real_,
          refit_date = dates[t_refit]
        )
        out_pos <- out_pos + 1
      }
      next
    }

    fc <- safe_forecast(fit, n_ahead = length(forecast_idx))

    if (is.null(fc)) {
      for (j in forecast_idx) {
        out_list[[out_pos]] <- data.frame(
          date = dates[j],
          mu = NA_real_,
          sigma = NA_real_,
          sigma2 = NA_real_,
          resid = NA_real_,
          z = NA_real_,
          refit_date = dates[t_refit]
        )
        out_pos <- out_pos + 1
      }
      next
    }

    mu_fc <- as.numeric(fitted(fc))
    sigma_fc <- as.numeric(sigma(fc))

    for (h in seq_along(forecast_idx)) {
      j <- forecast_idx[h]
      r_next <- ret[j]
      mu_h <- mu_fc[h]
      sigma_h <- sigma_fc[h]
      sigma2_h <- sigma_h^2

      z_next <- if (is.na(sigma_h) || sigma_h <= 0) {
        NA_real_
      } else {
        (r_next - mu_h) / sigma_h
      }

      out_list[[out_pos]] <- data.frame(
        date = dates[j],
        mu = mu_h,
        sigma = sigma_h,
        sigma2 = sigma2_h,
        resid = r_next - mu_h,
        z = z_next,
        refit_date = dates[t_refit]
      )
      out_pos <- out_pos + 1
    }
  }

  out <- do.call(rbind, out_list)

  # sort and return
  out <- out[order(out$date), ]
  rownames(out) <- NULL
  out
}

run_all_garch_monthly <- function(data, date_col = "date", init_window = 504) {
  dates <- data[[date_col]]
  ret_names <- setdiff(names(data), date_col)

  out <- lapply(ret_names, function(nm) {
    fit_garch_expanding_monthly(
      ret = data[[nm]],
      dates = dates,
      init_window = init_window
    )
  })

  names(out) <- ret_names
  out
}

# Too long to run what is below, so will run test case first
# garch_results <- run_all_garch_daily(managed_portfolios)

system.time({
  test_result <- fit_garch_expanding_monthly(
    ret = managed_portfolios[[2]][1:1100],
    dates = managed_portfolios$date[1:1100],
    init_window = 504
  )
})

system.time({
  garch_results_test <- run_all_garch_monthly(
    managed_portfolios[1:1100, ],
    date_col = "date",
    init_window = 504
  )
})

# Checking if there are any nas  and if no of forecasts is 596
sapply(garch_results_test, function(x) sum(is.na(x$sigma)))
sapply(garch_results_test, nrow)

# Constructing standardised residual matrix Z and Diagonal volatility matrix D, calculated
# in GARCH functions already
build_garch_blocks <- function(garch_results_list) {
  # Extracting reference dates from the first factors
  ref_dates <- garch_results_list[[1]]$date
  
  # Checking if all factors have identical dates
  same_dates <- all(sapply(garch_results_list, function(x) identical(x$date, ref_dates)))
  if (!same_dates) {
    stop("Dates are not aligned across assets.")
  }
  
  # Loop over each factor, extract standardised residuals, then put in a matrix
  Z <- sapply(garch_results_list, function(x) x$z)
  SIGMA <- sapply(garch_results_list, function(x) x$sigma)
  MU <- sapply(garch_results_list, function(x) x$mu)
  RESID <- sapply(garch_results_list, function(x) x$resid)
  
  Z <- as.matrix(Z)
  SIGMA <- as.matrix(SIGMA)
  MU <- as.matrix(MU)
  RESID <- as.matrix(RESID)
  
  # Setting column names
  colnames(Z) <- names(garch_results_list)
  colnames(SIGMA) <- names(garch_results_list)
  colnames(MU) <- names(garch_results_list)
  colnames(RESID) <- names(garch_results_list)
  
  # Setting each row to a trading day
  rownames(Z) <- as.character(ref_dates)
  rownames(SIGMA) <- as.character(ref_dates)
  rownames(MU) <- as.character(ref_dates)
  rownames(RESID) <- as.character(ref_dates)
  
  list(
    dates = ref_dates,
    Z = Z,
    SIGMA = SIGMA,
    MU = MU,
    RESID = RESID
  )
}

garch_blocks <- build_garch_blocks(garch_results_test)

Z_block <- garch_blocks$Z
SIGMA_block <- garch_blocks$SIGMA
MU_block <- garch_blocks$MU
RESID_block <- garch_blocks$RESID
dates_block <- garch_blocks$dates

# Dimension check
dim(Z_block)
dim(SIGMA_block)
dim(MU_block)
dim(RESID_block)

identical(dim(Z_block), dim(SIGMA_block))
identical(dim(Z_block), dim(MU_block))
identical(dim(Z_block), dim(RESID_block))

head(rownames(Z_block))
head(rownames(SIGMA_block))
head(colnames(Z_block))

# Missing value check
sum(is.na(Z_block))
sum(is.na(SIGMA_block))
sum(is.na(MU_block))
sum(is.na(RESID_block))

sum(!is.finite(Z_block))
sum(!is.finite(SIGMA_block))
sum(!is.finite(MU_block))
sum(!is.finite(RESID_block))

# Sigma positivity check
summary(as.vector(SIGMA_block))
sum(SIGMA_block <= 0)

# Identitiy check
max(abs(Z_block - (RESID_block / SIGMA_block)), na.rm = TRUE)

# Standardization checks (check if mean and sd are approx 0 and 1)
print(mean(as.vector(Z_block), na.rm = TRUE))
print(sd(as.vector(Z_block), na.rm = TRUE))
print(summary(colMeans(Z_block, na.rm = TRUE)))
print(summary(apply(Z_block, 2, sd, na.rm = TRUE)))

# Extreme values
print(max(abs(Z_block), na.rm = TRUE))
print(quantile(abs(Z_block), probs = c(0.90, 0.95, 0.99, 0.999), na.rm = TRUE))

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
