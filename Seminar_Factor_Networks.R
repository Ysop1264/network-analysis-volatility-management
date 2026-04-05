#### Uncomment if packages not installed
#install.packages("tidyverse")
#install.packages("tidyfinance")
#install.packages("scales")
#install.packages("frenchdata")
#install.packages("dplyr")
#install.packages("moments")
#install.packages("sandwich")
#install.packages("rlang")
#install.packages("lmtest")
#install.packages("lubridate")
#install.packages("nlshrink")
#install.packages("rugarch")
#install.packages("xts")
#install.packages("zoo")
#install.packages("igraph")
# install.packages("PeerPerformance")
####

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
library(igraph)
library(PeerPerformance)

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
    across(-c(date, mkt_excess, smb, hml, rmw, cma, risk_free, mom), ~ . - risk_free)
  )

managed_portfolios <- factors_joined_excess |> select(-risk_free)

print(head(managed_portfolios))

# =================================================
# Functions
# Functions
# =================================================

#' Compute the MVE static weights
#' 
#' 
#' @param mu Vector of means
#' @param sigma Variance-covariance matrix of returns
#' @return Solved weights b
compute_MVE_weights <- function(mu, sigma) {
  b_raw <- as.vector(solve(sigma, mu))
  b <- b_raw / sum(b_raw)
  return(b)
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

#'Creates the summary statistics table of managed portfolios returns.
#'The means, standard deviations are annualized and in percentages.
#'The min and max are in percentages.
#'The Sharpe Ratio is annualized.
#'
#'@param factor_returns data frame containing the returns of managed portfolios
#'@return returns the data frame with summary statistics
#'
summary_stats_table <- function(factor_returns){
  if(class(factor_returns[[1]]) == "Date"){
    factor_returns <- factor_returns[,-1, drop = FALSE]
  }
  
  means <- colMeans(factor_returns) * 100 * 252
  sd <- sapply(factor_returns, sd) * 100 * sqrt(252)
  kurtosis <- sapply(factor_returns, kurtosis)
  skewness <- sapply(factor_returns, skewness)
  min <- sapply(factor_returns, min) * 100
  max <- sapply(factor_returns, max) * 100
  Sharpe <- means/sd 
  
  df <- data.frame(
    Mean = round(means,4),
    SD = round(sd, 4),
    Kurtosis = round(kurtosis, 4),
    Skewness = round(skewness, 4),
    Min = round(min, 4),
    Max = round(max, 4),
    SR = round(Sharpe,4)
    
  )
  return(df)
}

#'Creates the short summary statistics table of managed portfolios returns.
#'The means, standard deviations are annualized and in percentages.
#'The min and max are in percentages.
#'The Sharpe Ratio is annualized.
#'
#'@param factor_returns data frame containing the returns of managed portfolios
#'@return returns the data frame with summary statistics
#'
summary_stats_short_table <- function(factor_returns){
  return(t(summary_stats_table(factor_returns)))
}

#' Graphs the standard line plot of annualized returns
#'
#'@param returns_df data frame containing date as first column, and returns as second column
#'
graph_annualized_returns <- function(returns_df){
  
  # Changing how often there are ticks on xaxis and y axis
  axis <- par(lab = c(20, 8, 5))
  # Plotting the returns 
  plot(main = "", xlab = "Date", ylab = "Annaulized returns",
       y = returns_df[[2]] * 252 , x = as.Date(returns_df[[1]]) , type = "l" , ylim = c(-40,40))
  
  # Creates ticks with rounded values of the annualized returns 
  y_ticks <- pretty(returns_df[[2]] * 252)
  abline(h = y_ticks, col = "grey85", lty = 1)
  lines(y = returns_df[[2]] * 252 , x = as.Date(returns_df[[1]]))
  
  
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
  init_window = 504,
  garch_order = c(1, 1),
  arma_order = c(0, 0),
  distribution = "norm",
  scale_ret = TRUE,
  series_name = NA_character_
) {
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

  safe_fit <- function(x, series_name = NA, refit_date = NA) {
    tryCatch(
      suppressWarnings(
        ugarchfit(
          spec = spec,
          data = x,
          solver = "hybrid",
          solver.control = list(trace = 0)
        )
      ),
      error = function(e) {
        message("FIT FAILED | series=", series_name,
                " | refit_date=", refit_date,
                " | msg=", conditionMessage(e))
        NULL
      }
    )
  }

  safe_forecast <- function(fit, n_ahead, series_name = NA, refit_date = NA) {
    tryCatch(
      ugarchforecast(fit, n.ahead = n_ahead),
      error = function(e) {
        message("FORECAST FAILED | series=", series_name,
                " | refit_date=", refit_date,
                " | msg=", conditionMessage(e))
        NULL
      }
    )
  }

  idx_all <- seq_len(n)
  month_id <- format(dates, "%Y-%m")

  month_end_idx <- tapply(idx_all, month_id, max)
  month_end_idx <- as.integer(month_end_idx)

  refit_idx <- month_end_idx[month_end_idx >= init_window & month_end_idx < n]

  if (length(refit_idx) == 0) {
    stop("No eligible monthly refit dates found after init_window.")
  }

  out_list <- vector("list", n - init_window)
  out_pos <- 1
  
  # Fit only once 
  t_refit <- init_window
  insample_ret <- ret[1:t_refit]
  
  fit <- safe_fit(insample_ret, series_name, dates[t_refit])
  
  if (is.null(fit)) {
    stop("Initial GARCH fit failed — cannot freeze parameters.")
  }
  
  # extract fixed parameters
  fixed_pars <- coef(fit)
  
  # create fixed spec
  spec_fixed <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = garch_order),
    mean.model = list(armaOrder = arma_order, include.mean = TRUE),
    distribution.model = distribution,
    fixed.pars = as.list(fixed_pars)
  )
  
    # forecast entire remaining sample with fixed parameters
  fc <- ugarchforecast(
    spec_fixed,
    data = ret,
    out.sample = n - t_refit,
    n.ahead = 1,
    n.roll = n - t_refit - 1
  )

  mu_fc <- as.numeric(fitted(fc))
  sigma_fc <- as.numeric(sigma(fc))

  n_oos <- n - t_refit
  if (length(mu_fc) != n_oos || length(sigma_fc) != n_oos) {
    stop("Forecast length does not match out-of-sample length.")
  }

  out_list <- vector("list", n_oos)

  for (j in seq_len(n_oos)) {
    idx <- t_refit + j
    r_next <- ret[idx]

    sigma_h <- sigma_fc[j]
    mu_h <- mu_fc[j]

    z_next <- if (is.na(sigma_h) || sigma_h <= 0) {
      NA_real_
    } else {
      (r_next - mu_h) / sigma_h
    }

    out_list[[j]] <- data.frame(
      date = dates[idx],
      mu = mu_h,
      sigma = sigma_h,
      sigma2 = sigma_h^2,
      resid = r_next - mu_h,
      z = z_next,
      refit_date = dates[t_refit]
    )
  }
  
  out <- do.call(rbind, out_list)
  out <- out[order(out$date), ]
  rownames(out) <- NULL
  return(out)
  
  # comment out for now to freeze the estimates 2 years 
  # for (k in seq_along(refit_idx)) {
  #   t_refit <- refit_idx[k]
  #   insample_ret <- ret[1:t_refit]
  # 
  #   next_refit <- if (k < length(refit_idx)) refit_idx[k + 1] else n
  #   forecast_idx <- (t_refit + 1):next_refit
  # 
  #   if (length(forecast_idx) == 0) next
  # 
  #   if (sd(insample_ret, na.rm = TRUE) <= 1e-8) {
  #     for (j in forecast_idx) {
  #       out_list[[out_pos]] <- data.frame(
  #         date = dates[j],
  #         mu = NA_real_,
  #         sigma = NA_real_,
  #         sigma2 = NA_real_,
  #         resid = NA_real_,
  #         z = NA_real_,
  #         refit_date = dates[t_refit]
  #       )
  #       out_pos <- out_pos + 1
  #     }
  #     next
  #   }
  # 
  #   fit <- safe_fit(insample_ret, series_name, dates[t_refit])
  # 
  #   if (is.null(fit)) {
  #     for (j in forecast_idx) {
  #       out_list[[out_pos]] <- data.frame(
  #         date = dates[j],
  #         mu = NA_real_,
  #         sigma = NA_real_,
  #         sigma2 = NA_real_,
  #         resid = NA_real_,
  #         z = NA_real_,
  #         refit_date = dates[t_refit]
  #       )
  #       out_pos <- out_pos + 1
  #     }
  #     next
  #   }
  # 
  #   fc <- safe_forecast(fit, length(forecast_idx), series_name, dates[t_refit])
  # 
  #   if (!is.null(fit)) {
  #     message("PARAM NAMES | series=", series_name, " | ",
  #             paste(names(coef(fit)), collapse = ", "))
  #   }
  # 
  #   if (is.null(fc)) {
  #     for (j in forecast_idx) {
  #       out_list[[out_pos]] <- data.frame(
  #         date = dates[j],
  #         mu = NA_real_,
  #         sigma = NA_real_,
  #         sigma2 = NA_real_,
  #         resid = NA_real_,
  #         z = NA_real_,
  #         refit_date = dates[t_refit]
  #       )
  #       out_pos <- out_pos + 1
  #     }
  #     next
  #   }
  # 
  #   mu_fc <- as.numeric(fitted(fc))
  #   sigma_fc <- as.numeric(sigma(fc))
  # 
  #   for (h in seq_along(forecast_idx)) {
  #     j <- forecast_idx[h]
  #     r_next <- ret[j]
  #     mu_h <- mu_fc[h]
  #     sigma_h <- sigma_fc[h]
  #     sigma2_h <- sigma_h^2
  # 
  #     z_next <- if (is.na(sigma_h) || sigma_h <= 0) {
  #       NA_real_
  #     } else {
  #       (r_next - mu_h) / sigma_h
  #     }
  # 
  #     out_list[[out_pos]] <- data.frame(
  #       date = dates[j],
  #       mu = mu_h,
  #       sigma = sigma_h,
  #       sigma2 = sigma2_h,
  #       resid = r_next - mu_h,
  #       z = z_next,
  #       refit_date = dates[t_refit]
  #     )
  #     out_pos <- out_pos + 1
  #   }
  # }

  # out <- do.call(rbind, out_list)
  # out <- out[order(out$date), ]
  # rownames(out) <- NULL
  # out
}

run_all_garch_monthly <- function(data, date_col = "date", init_window = 504) {
  dates <- data[[date_col]]
  ret_names <- setdiff(names(data), date_col)

  out <- lapply(ret_names, function(nm) {
    fit_garch_expanding_monthly(
      ret = data[[nm]],
      dates = dates,
      init_window = init_window,
      series_name = nm,
      distribution = "std"
    )
  })

  names(out) <- ret_names
  out
}

# Too long to run what is below, so will run test case first
# garch_results <- run_all_garch_daily(managed_portfolios)

#system.time({
 # test_result <- fit_garch_expanding_monthly(
  #  ret = managed_portfolios[[2]][1:1100],
   # dates = managed_portfolios$date[1:1100],
    #init_window = 504,
   # series_name = names(managed_portfolios)[2]
  #)
#})

#sum(is.na(test_result$sigma))
#nrow(test_result)
#head(test_result)
#tail(test_result)

#system.time({
#  test_result <- fit_garch_expanding_monthly(
#    ret = managed_portfolios[[2]][1:1100],
#    dates = managed_portfolios$date[1:1100],
#    init_window = 504
#  )
#})

#system.time({
#  garch_results_test <- run_all_garch_monthly(
#    managed_portfolios[1:1100, ],
#    date_col = "date",
#    init_window = 504
#  )
#})

#na_summary <- sapply(garch_results_test, function(x) sum(is.na(x$sigma)))
#sort(na_summary[na_summary > 0], decreasing = TRUE)

#table(sapply(garch_results_test, nrow))

#z_summary <- sapply(garch_results_test, function(x) sd(x$z, na.rm = TRUE))
#summary(z_summary)

# sGARCH is underestimating volatility (has fat tails , sd of 1.18, 
# using Studnet t-distribution for a run)
#system.time({
#  garch_results_test_t <- run_all_garch_monthly(
#    managed_portfolios[1:1100, ],
#    date_col = "date",
#    init_window = 504
#  )
#})

#na_summary_t <- sapply(garch_results_test_t, function(x) sum(is.na(x$sigma)))
#sort(na_summary_t[na_summary_t > 0], decreasing = TRUE)

#z_summary_t <- sapply(garch_results_test_t, function(x) sd(x$z, na.rm = TRUE))
#summary(z_summary_t)

#table(sapply(garch_results_test_t, nrow))

#mean_z <- sapply(garch_results_test_t, function(x) mean(x$z, na.rm = TRUE))
#summary(mean_z)

# Constructing standardised residual matrix Z, calculated
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
  Z <- do.call(cbind, lapply(garch_results_list, `[[`, "z"))
  SIGMA <- do.call(cbind, lapply(garch_results_list, `[[`, "sigma"))
  MU <- do.call(cbind, lapply(garch_results_list, `[[`, "mu"))
  RESID <- do.call(cbind, lapply(garch_results_list, `[[`, "resid"))
  
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

# Keeping track of current month and the forecasting month
# Minimum estimation window of 1 trading year, or else DCC becomes too unstable
make_month_index <- function(dates, min_obs = 252, rolling_window = NULL) {
  dates <- as.Date(dates)
  n <- length(dates)

  if (n == 0) stop("dates is empty.")

  month_id <- format(dates, "%Y-%m")

  # Finding the last trading day of each month
  month_end_idx <- tapply(seq_len(n), month_id, max)
  refit_months <- names(month_end_idx)
  month_end_idx <- as.integer(month_end_idx)

  # the final month has no following month to forecast
  out <- vector("list", length(month_end_idx) - 1)
  out_pos <- 1

  for (k in seq_len(length(month_end_idx) - 1)) {
    t_end <- month_end_idx[k]
    next_end <- month_end_idx[k + 1]

    if (is.null(rolling_window)) {
      insample_idx <- seq_len(t_end)
    } else {
      start_idx <- max(1, t_end - rolling_window + 1)
      insample_idx <- start_idx:t_end
    }

    forecast_idx <- (t_end + 1):next_end

    if (length(insample_idx) < min_obs) next
    if (length(forecast_idx) == 0) next

    out[[out_pos]] <- list(
      refit_month = refit_months[k],
      refit_idx = t_end,
      refit_date = dates[t_end],
      insample_idx = insample_idx,
      forecast_idx = forecast_idx,
      forecast_dates = dates[forecast_idx]
    )
    out_pos <- out_pos + 1
  }

  out <- out[seq_len(out_pos - 1)]
  out
}

# Function to convert covariance matrix to correlation
cov_to_cor <- function(Sigma, eps = 1e-10) {
  d <- sqrt(pmax(diag(Sigma), eps))
  Dinv <- diag(1 / d, nrow = length(d))
  Dinv %*% Sigma %*% Dinv
}

# Forcing matrix to be positive semi-definite
# Done to allow inversion, valid likelihood and prevent DCC explosion
# Basically trying to find the closeest valid covariance matrix
make_psd <- function(M, eps = 1e-8) {
  M <- (M + t(M)) / 2
  
  # Eigenvector decomposition
  ee <- eigen(M, symmetric = TRUE)

  # Pushing any negative eigenvalues to small positive number
  vals <- pmax(ee$values, eps)
  out <- ee$vectors %*% diag(vals, nrow = length(vals)) %*% t(ee$vectors)
  out <- (out + t(out)) / 2
  out
}

# Constructs robust correlation target matrix S
safe_cor <- function(Z) {
  S <- suppressWarnings(cor(Z, use = "pairwise.complete.obs"))
  S[!is.finite(S)] <- 0
  diag(S) <- 1
  S <- make_psd(S)
  S <- cov_to_cor(S)
  S
}

# Nonlinear-shrinkage long-run target for DCC
# Input: in-sample standardized residuals Z_insample
# Output: correlation target matrix S_target
build_nlshrink_target <- function(Z_insample, eps = 1e-8) {
  Z_insample <- as.matrix(Z_insample)

  # DCC target must be built on complete rows only
  Z_use <- Z_insample[complete.cases(Z_insample), , drop = FALSE]

  if (nrow(Z_use) < 20) {
    stop("Too few complete observations for nonlinear shrinkage.")
  }

  # Standardized residuals should already be centered approximately at zero,
  # but demeaning helps the covariance estimator
  Z_use <- scale(Z_use, center = TRUE, scale = FALSE)

  # nonlinear shrinkage covariance estimate
  Sigma_nl <- nlshrink::nlshrink_cov(Z_use)

  # Safety cleanup
  Sigma_nl <- make_psd(Sigma_nl, eps = eps)

  # Convert shrunk covariance target into a correlation target for DCC
  S_target <- cov_to_cor(Sigma_nl, eps = eps)
  # Correlation might not be PSD, so re-done even though Sigma_nl is psd
  S_target <- make_psd(S_target, eps = eps)
  # Now, diagonals are no longer 1, so making it a corr matrix again
  S_target <- cov_to_cor(S_target, eps = eps)

  S_target
}

# Converts DCC matrix Q_t to correlation matrix R_t
normalize_Q_to_R <- function(Q, eps = 1e-10) {
  qd <- sqrt(pmax(diag(Q), eps))
  Dinv <- diag(1 / qd, nrow = length(qd))
  R <- Dinv %*% Q %*% Dinv
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R
}

# Function to run the DCC recursion over the time series
dcc_filter <- function(Z, a, b, S) {
  Z <- as.matrix(Z)
  Tn <- nrow(Z) # number of time obs
  N <- ncol(Z) # number of factors and managed portfolios

  # Store the entire time series of intermediate matrix and correlation matrix
  Q_list <- vector("list", Tn)
  R_list <- vector("list", Tn)

  # Initialising recursion at the long run target
  Q_prev <- S

  for (t in seq_len(Tn)) {
    # for DCC recursion at time t, use lagged standardised residual t-1 (at first, just use 0)
    zlag <- if (t == 1) rep(0, N) else Z[t - 1, ]
    zlag[!is.finite(zlag)] <- 0

    # DCC Recursion (tcrossprod is tranpose product)
    Q_t <- (1 - a - b) * S + a * tcrossprod(zlag) + b * Q_prev
    # Enforce symmetry to make it diagonalisable
    Q_t <- (Q_t + t(Q_t)) / 2
    # Convert covariance object to correlation matrix
    R_t <- normalize_Q_to_R(Q_t)

    Q_list[[t]] <- Q_t
    R_list[[t]] <- R_t
    Q_prev <- Q_t
  }

  list(Q = Q_list, R = R_list, Q_T = Q_prev)
}

# Function computes the Gaussian Negative Likelihood for the DCC Model
# for given parameters (a, b)
dcc_negloglik <- function(par, Z, S, penalty = 1e12, ridge = 1e-8) {
  a <- par[1]
  b <- par[2]

  # Constraint check to ensure feasibility
  if (!is.finite(a) || !is.finite(b) || a < 0 || b < 0 || (a + b) >= 0.999) {
    return(penalty)
  }

  Z <- as.matrix(Z)
  S <- as.matrix(S)

  Tn <- nrow(Z)
  N <- ncol(Z)

  Q_prev <- S
  ll <- 0

  for (t in seq_len(Tn)) {
    zlag <- if (t == 1) rep(0, N) else Z[t - 1, ]
    zlag[!is.finite(zlag)] <- 0

    Q_t <- (1 - a - b) * S + a * tcrossprod(zlag) + b * Q_prev
    Q_t <- (Q_t + t(Q_t)) / 2

    Rt <- normalize_Q_to_R(Q_t, eps = ridge)
    if (any(!is.finite(Rt))) {
      return(penalty)
    }

    Rt <- make_psd(Rt, eps = ridge)
    Rt <- (Rt + t(Rt)) / 2
    diag(Rt) <- 1

    # Ridge regularisation to avoid near singular correlation matrices
    Rt_reg <- Rt + diag(ridge, N)

    # Cholesky Decomposition for faster output
    U <- tryCatch(chol(Rt_reg), error = function(e) NULL)
    if (is.null(U)) {
      return(penalty)
    }

    logdet <- 2 * sum(log(diag(U)))
    if (!is.finite(logdet)) {
      return(penalty)
    }

    zt <- as.numeric(Z[t, ])
    if (any(!is.finite(zt))) {
      return(penalty)
    }

    # backsolve solves the triangular factor in the Cholesky without forming the inverse directly
    y <- tryCatch(backsolve(U, zt, transpose = TRUE), error = function(e) NULL)
    if (is.null(y) || any(!is.finite(y))) {
      return(penalty)
    }

    x <- tryCatch(backsolve(U, y), error = function(e) NULL)
    if (is.null(x) || any(!is.finite(x))) {
      return(penalty)
    }

    quad <- sum(zt * x)
    if (!is.finite(quad)) {
      return(penalty)
    }

    ll <- ll + logdet + quad
    Q_prev <- Q_t
  }

  0.5 * ll
}

# Wrapper function that does the DCC estimation
estimate_dcc <- function(Z_insample, S_target = NULL, start_par = c(0.03, 0.95), # typical small a, large b
                         ridge = 1e-8, stationarity_eps = 1e-4) {
  Z_insample <- as.matrix(Z_insample)
  # Missing values distort recursion and the likelihood
  Z_use <- Z_insample[complete.cases(Z_insample), , drop = FALSE]

  if (nrow(Z_use) < 50) {
    stop("Too few complete observations to estimate DCC.")
  }

  # Falls back to robust correlation target if no shrinkage target is supplied
  if (is.null(S_target)) {
    S_target <- safe_cor(Z_use)
  } else {
    S_target <- as.matrix(S_target)
    S_target <- (S_target + t(S_target)) / 2
    S_target <- make_psd(S_target)
    S_target <- cov_to_cor(S_target)
    S_target <- make_psd(S_target)
    diag(S_target) <- 1
  }

  start_a <- max(1e-6, start_par[1])
  start_b <- max(1e-6, start_par[2])

  if ((start_a + start_b) >= (1 - stationarity_eps)) {
    scale_factor <- (1 - stationarity_eps) / (start_a + start_b)
    start_a <- start_a * scale_factor
    start_b <- start_b * scale_factor
  }

  # Uses box constraints to keep parameters positive and from wandering too far
  opt <- tryCatch(
    optim(
      par = c(start_a, start_b),
      fn = function(par) dcc_negloglik(par, Z_use, S_target, ridge = ridge),
      method = "L-BFGS-B",
      lower = c(1e-6, 1e-6),
      upper = c(0.5, 0.999)
    ),
    error = function(e) NULL
  )

  if (is.null(opt) ||
      !is.list(opt) ||
      is.null(opt$convergence) ||
      opt$convergence != 0 ||
      is.null(opt$value) ||
      !is.finite(opt$value)) {
    stop("DCC optimization failed to converge.")
  }

  a_hat <- opt$par[1]
  b_hat <- opt$par[2]

  if ((a_hat + b_hat) >= (1 - stationarity_eps)) {
    scale_factor <- (1 - stationarity_eps) / (a_hat + b_hat)
    a_hat <- a_hat * scale_factor
    b_hat <- b_hat * scale_factor
  }

  # To obtain full in-sample Q_t, full in-sample R_t, and terminal Q_T
  filt <- dcc_filter(Z_use, a_hat, b_hat, S_target)

  list(
    a = a_hat,
    b = b_hat,
    S = S_target,
    Q_T = filt$Q_T,
    Q = filt$Q,
    R = filt$R,
    nll = opt$value,
    convergence = opt$convergence
  )
}

# Function that uses DCC estimation to forecast multiple steps ahead correlations
forecast_dcc_correlations <- function(Q_T, S, a, b, h ) {
  Q_T <- as.matrix(Q_T)
  S <- as.matrix(S)

  Q_T <- (Q_T + t(Q_T)) / 2
  S <- (S + t(S)) / 2

  Q_T <- make_psd(Q_T)
  S <- make_psd(S)
  S <- cov_to_cor(S)
  S <- make_psd(S)
  diag(S) <- 1

  phi <- a + b
  if (!is.finite(phi) || phi < 0 || phi >= 1) {
    stop("Invalid DCC persistence parameter: a + b must lie in [0, 1).")
  }

  Q_fc <- vector("list", h)
  R_fc <- vector("list", h)

  Q_prev <- Q_T

  for (k in seq_len(h)) {
    Q_next <- S + phi * (Q_prev - S)
    Q_next <- (Q_next + t(Q_next)) / 2
    Q_next <- make_psd(Q_next)

    R_next <- normalize_Q_to_R(Q_next)
    R_next <- make_psd(R_next)
    R_next <- cov_to_cor(R_next)
    diag(R_next) <- 1

    Q_fc[[k]] <- Q_next
    R_fc[[k]] <- R_next
    Q_prev <- Q_next
  }

  list(Q = Q_fc, R = R_fc)
}

# Rebuilding conditional covariance matrices using future volatility
# forecasts from univariate GARCH and future correlation forecasts from DCC
build_cov_from_sigma_and_R <- function(SIGMA_future, R_list, dates_future = NULL) {
  SIGMA_future <- as.matrix(SIGMA_future) # Matrix of daily forecats over the future month
  h <- nrow(SIGMA_future)
  N <- ncol(SIGMA_future)

  if (length(R_list) != h) {
    stop("length(R_list) must equal nrow(SIGMA_future).")
  }

  H_list <- vector("list", h)

  for (t in seq_len(h)) {
    sig_t <- SIGMA_future[t, ]
    if (any(!is.finite(sig_t)) || any(sig_t <= 0)) {
      H_list[[t]] <- matrix(NA_real_, N, N)
      next
    }

    # Main covariance reconstruction formula
    D_t <- diag(sig_t, nrow = N)
    R_t <- R_list[[t]]
    H_t <- D_t %*% R_t %*% D_t
    H_t <- (H_t + t(H_t)) / 2
    H_list[[t]] <- H_t
  }

  if (!is.null(dates_future)) {
    names(H_list) <- as.character(as.Date(dates_future))
  }

  # Returns a list of daily covariance forecast matrices
  H_list
}

# Convert daily covariance forecasts to monthly
aggregate_monthly_cov <- function(H_list) {
  valid <- Filter(function(x) all(is.finite(x)), H_list)
  if (length(valid) == 0) return(NULL)

  Reduce(`+`, valid)
}

# Grouping helpers (we finally decided on biyearly)
quarter_id <- function(x) {
  x <- as.Date(x)
  yr <- format(x, "%Y")
  qtr <- ((as.integer(format(x, "%m")) - 1) %/% 3) + 1
  paste0(yr, "-Q", qtr)
}

year_id <- function(x) {
  x <- as.Date(x)
  format(x, "%Y")
}

biyear_id <- function(x) {
  x <- as.Date(x)
  yr <- as.integer(format(x, "%Y"))
  start_yr <- ifelse(yr %% 2 == 0, yr - 1, yr)
  paste0(start_yr, "-", start_yr + 1)
}

# Master function
# loops over month-end refit dates and does the full monthly DCC forecasting exercise
# estimate DCC on expanding in-sample data
# forecast next-month daily correlations
# combine with next-month daily volatility forecasts
# aggregate into monthly covariance forecasts
run_dcc_monthly <- function(Z_block, SIGMA_block, dates_block,
                            min_corr_window = 252,
                            rolling_window = 504,
                            S_builder = NULL,
                            S_update_freq = c("monthly", "quarterly", "yearly", "biyearly"),
                            trace = TRUE) {
  frozen_ab <- NULL
  frozen_S <- NULL
  Z_block <- as.matrix(Z_block)
  SIGMA_block <- as.matrix(SIGMA_block)
  dates_block <- as.Date(dates_block)
  S_update_freq <- match.arg(S_update_freq)

  if (!identical(dim(Z_block), dim(SIGMA_block))) {
    stop("Z_block and SIGMA_block must have identical dimensions.")
  }

  if (nrow(Z_block) != length(dates_block)) {
    stop("dates_block length must equal nrow(Z_block).")
  }

  month_map <- make_month_index(
    dates = dates_block,
    min_obs = min_corr_window,
    rolling_window = rolling_window
  )

  if (length(month_map) == 0) {
    stop("No eligible monthly DCC forecasting blocks found.")
  }

  out <- vector("list", length(month_map))

  last_S_target <- NULL
  last_S_period <- NULL

  for (i in seq_along(month_map)) {
    blk <- month_map[[i]]

    if (trace) {
      message(
        "Estimating DCC for refit date ",
        as.character(blk$refit_date),
        " | insample n = ",
        length(blk$insample_idx),
        " | forecast h = ",
        length(blk$forecast_idx)
      )
    }

    Z_insample <- Z_block[blk$insample_idx, , drop = FALSE]
    SIGMA_future <- SIGMA_block[blk$forecast_idx, , drop = FALSE]

    current_S_period <- switch(
      S_update_freq,
      monthly   = format(blk$refit_date, "%Y-%m"),
      quarterly = quarter_id(blk$refit_date),
      yearly    = year_id(blk$refit_date),
      biyearly  = biyear_id(blk$refit_date)
    )

    reused_S <- FALSE

  #   if (is.null(S_builder)) {
  #     S_target <- safe_cor(Z_insample)
  #     } else {
  #   if (is.null(last_S_target) || !identical(current_S_period, last_S_period)) {
  #   if (trace) {
  #     message("Updating S_target via ", S_update_freq,
  #             " NL shrinkage at ", as.character(blk$refit_date))
  #   }
  #   S_target <- S_builder(Z_insample)
  #   last_S_target <- S_target
  #   last_S_period <- current_S_period
  # } else {
  #   S_target <- last_S_target
  #   reused_S <- TRUE
  # }

  if (!is.null(frozen_S)) {
  S_target <- frozen_S
  reused_S <- TRUE
  } else if (is.null(S_builder)) {
    S_target <- safe_cor(Z_insample)
  } else {
    if (is.null(last_S_target) || !identical(current_S_period, last_S_period)) {
      if (trace) {
        message("Updating S_target via ", S_update_freq,
              " NL shrinkage at ", as.character(blk$refit_date))
      }
    S_target <- S_builder(Z_insample)
    last_S_target <- S_target
    last_S_period <- current_S_period
  } else {
    S_target <- last_S_target
    reused_S <- TRUE
  }
}

    cat("\n============================\n")
    cat("Refit date:", as.character(blk$refit_date), "\n")
    cat("dim(Z_insample):", paste(dim(Z_insample), collapse = " x "), "\n")
    cat("Any non-finite in Z_insample:", any(!is.finite(Z_insample)), "\n")
    cat("Complete rows:", sum(complete.cases(Z_insample)), "\n")

    z_sds <- apply(Z_insample, 2, sd, na.rm = TRUE)
    cat("Min column sd:", min(z_sds, na.rm = TRUE), "\n")
    cat("Any zero/NA column sd:",
        any(!is.finite(z_sds) | z_sds < 1e-8), "\n")

    cat("S update period:", current_S_period, "\n")
    cat("Reused previous S_target:", reused_S, "\n")

    S_dbg <- S_target
    cat("Any non-finite in S:", any(!is.finite(S_dbg)), "\n")
    cat("Symmetry check:", max(abs(S_dbg - t(S_dbg))), "\n")
    cat("Min eigenvalue of S:",
        min(eigen(S_dbg, symmetric = TRUE, only.values = TRUE)$values), "\n")
    cat("============================\n")
    
    # Comment out to try freezing estimates 
    # dcc_fit <- tryCatch(
    #   estimate_dcc(Z_insample, S_target),
    #   error = function(e) {
    #     cat("DCC failed at", as.character(blk$refit_date), ":", e$message, "\n")
    #     return(NULL)
    #   }
    # )
    if (is.null(frozen_ab)) {
  dcc_fit <- tryCatch(
    estimate_dcc(Z_insample, S_target),
    error = function(e) {
      cat("DCC failed at", as.character(blk$refit_date), ":", e$message, "\n")
      return(NULL)
    }
  )

  if (!is.null(dcc_fit)) {
    frozen_ab <- c(dcc_fit$a, dcc_fit$b)
    frozen_S  <- dcc_fit$S
  }

} else {
  dcc_fit <- list(
    a = frozen_ab[1],
    b = frozen_ab[2],
    S = frozen_S,
    Q_T = dcc_filter(Z_insample, frozen_ab[1], frozen_ab[2], frozen_S)$Q_T
  )
}

if (is.null(dcc_fit)) {
  out[[i]] <- list(
    refit_date = blk$refit_date,
    forecast_dates = blk$forecast_dates,
    a = NA_real_,
    b = NA_real_,
    R_forecasts = NULL,
    Q_forecasts = NULL,
    H_forecasts = NULL,
    H_month = NULL
  )
  next
}

 dcc_fc <- forecast_dcc_correlations(
    Q_T = dcc_fit$Q_T,
    S = dcc_fit$S,
    a = dcc_fit$a,
    b = dcc_fit$b,
    h = length(blk$forecast_idx)
  )

  H_fc <- build_cov_from_sigma_and_R(
    SIGMA_future = SIGMA_future,
    R_list = dcc_fc$R,
    dates_future = blk$forecast_dates
  )

  H_month <- (1/length(blk$forecast_idx)) * aggregate_monthly_cov(H_fc) 

  out[[i]] <- list(
    refit_month = blk$refit_month,
    refit_date = blk$refit_date,
    insample_idx = blk$insample_idx,
    forecast_idx = blk$forecast_idx,
    forecast_dates = blk$forecast_dates,
    S_period = current_S_period,
    a = dcc_fit$a,
    b = dcc_fit$b,
    S = dcc_fit$S,
    Q_T = dcc_fit$Q_T,
    R_forecasts = dcc_fc$R,
    Q_forecasts = dcc_fc$Q,
    H_forecasts = H_fc,
    H_month = H_month
  )
  }

  names(out) <- sapply(out, function(x) as.character(x$refit_date))
  out
}

# test_data_all_medium <- managed_portfolios %>%
#   slice(1:1800)

# cat("Rows in medium test:", nrow(test_data_all_medium), "\n")
# cat("Number of assets/factors:", ncol(test_data_all_medium) - 1, "\n")

# time_garch_all_medium <- system.time({
#   garch_all_medium <- run_all_garch_monthly(
#     test_data_all_medium,
#     date_col = "date",
#     init_window = 504
#   )
# })

# print(time_garch_all_medium)

# blocks_all_medium <- build_garch_blocks(garch_all_medium)

# time_dcc_all_medium <- system.time({
#   dcc_all_medium <- run_dcc_monthly(
#     Z_block = blocks_all_medium$Z,
#     SIGMA_block = blocks_all_medium$SIGMA,
#     dates_block = blocks_all_medium$dates,
#     min_corr_window = 252,
#     rolling_window = 504,
#     S_builder = build_nlshrink_target,
#     S_update_freq = "monthly",
#     trace = TRUE
#   )
# })

# print(time_dcc_all_medium)

all_data_full <- managed_portfolios

cat("Rows in full sample:", nrow(all_data_full), "\n")
cat("Number of assets/factors:", ncol(all_data_full) - 1, "\n")

time_garch_all_full <- system.time({
  garch_all_full <- run_all_garch_monthly(
    all_data_full,
    date_col = "date",
    init_window = 504
  )
})

print(time_garch_all_full)

blocks_all_full <- build_garch_blocks(garch_all_full)

time_dcc_all_full <- system.time({
  dcc_all_full <- run_dcc_monthly(
    Z_block = blocks_all_full$Z,
    SIGMA_block = blocks_all_full$SIGMA,
    dates_block = blocks_all_full$dates,
    min_corr_window = 252,
    rolling_window = 504,
    S_builder = build_nlshrink_target,
    S_update_freq = "monthly",
    trace = TRUE
  )
})

print(time_dcc_all_full)

# Subset test for time check (fama french, momentum and random)
# test_data_cmp <- managed_portfolios %>%
#   select(
#     date,
#     mkt_excess,
#     smb,
#     hml,
#     rmw,
#     cma,
#     mom,
#     agric,     
#     food       
#   ) %>%
#   slice(1:1100)

# time_garch_cmp <- system.time({
#   garch_cmp <- run_all_garch_monthly(
#     test_data_cmp,
#     date_col = "date",
#     init_window = 504
#   )
# })

# print(time_garch_cmp)

# time_blocks_cmp <- system.time({
#   blocks_cmp <- build_garch_blocks(garch_cmp)
# })

# print(time_blocks_cmp)

# DONT RUN (NOT NEEDED)
# time_dcc_quarterly <- system.time({
#  dcc_quarterly <- run_dcc_monthly(
#    Z_block = blocks_cmp$Z,
#    SIGMA_block = blocks_cmp$SIGMA,
#    dates_block = as.Date(rownames(blocks_cmp$Z)),
#    min_corr_window = 252,
#    rolling_window = 504,
#    S_builder = build_nlshrink_target,
#    S_update_freq = "quarterly",
#    trace = TRUE
#  )
#})

#print(time_dcc_quarterly)

#time_dcc_yearly <- system.time({
#  dcc_yearly <- run_dcc_monthly(
#    Z_block = blocks_cmp$Z,
#    SIGMA_block = blocks_cmp$SIGMA,
#    dates_block = as.Date(rownames(blocks_cmp$Z)),
#    min_corr_window = 252,
#    rolling_window = 504,
#    S_builder = build_nlshrink_target,
#    S_update_freq = "yearly",
#    trace = TRUE
#  )
#})

#print(time_dcc_yearly)

# time_dcc_biyearly <- system.time({
#   dcc_biyearly <- run_dcc_monthly(
#     Z_block = blocks_cmp$Z,
#     SIGMA_block = blocks_cmp$SIGMA,
#     dates_block = as.Date(rownames(blocks_cmp$Z)),
#     min_corr_window = 252,
#     rolling_window = 504,
#     S_builder = build_nlshrink_target,
#     S_update_freq = "biyearly",
#     trace = TRUE
#   )
# })

# print(time_dcc_biyearly)


# get_dcc_param_summary <- function(dcc_obj) {
#   a_vals <- sapply(dcc_obj, function(x) x$a)
#   b_vals <- sapply(dcc_obj, function(x) x$b)

#   c(
#     mean_a = mean(a_vals, na.rm = TRUE),
#     mean_b = mean(b_vals, na.rm = TRUE),
#     mean_ab = mean(a_vals + b_vals, na.rm = TRUE),
#     median_a = median(a_vals, na.rm = TRUE),
#     median_b = median(b_vals, na.rm = TRUE),
#     median_ab = median(a_vals + b_vals, na.rm = TRUE)
#   )
# }

# cmp_table <- rbind(
#   quarterly = c(
#     elapsed = time_dcc_quarterly["elapsed"],
#     get_dcc_param_summary(dcc_quarterly)
#   ),
#   yearly = c(
#     elapsed = time_dcc_yearly["elapsed"],
#     get_dcc_param_summary(dcc_yearly)
#   ),
#   biyearly = c(
#     elapsed = time_dcc_biyearly["elapsed"],
#     get_dcc_param_summary(dcc_biyearly)
#   )
# )

# print(round(cmp_table, 4))

# average_H_month_matrix <- function(dcc_res) {
#   H_all <- lapply(dcc_res, `[[`, "H_month")
#   n <- nrow(H_all[[1]])
#   H_avg <- matrix(0, n, n)

#   for (k in seq_along(H_all)) {
#     H_avg <- H_avg + H_all[[k]]
#   }

#   H_avg / length(H_all)
# }

# H_q <- average_H_month_matrix(dcc_quarterly)
# H_y <- average_H_month_matrix(dcc_yearly)
# H_b <- average_H_month_matrix(dcc_biyearly)

# frobenius_dist <- function(A, B) {
#   sqrt(sum((A - B)^2))
# }

# comparison_H <- c(
#   q_vs_y = frobenius_dist(H_q, H_y),
#   q_vs_b = frobenius_dist(H_q, H_b),
#   y_vs_b = frobenius_dist(H_y, H_b)
# )

# print(comparison_H)

# monthwise_H_diff <- function(dcc1, dcc2) {
#   stopifnot(length(dcc1) == length(dcc2))

#   sapply(seq_along(dcc1), function(i) {
#     H1 <- dcc1[[i]]$H_month
#     H2 <- dcc2[[i]]$H_month
#     sqrt(sum((H1 - H2)^2))
#   })
# }

# diff_q_y <- monthwise_H_diff(dcc_quarterly, dcc_yearly)
# diff_q_b <- monthwise_H_diff(dcc_quarterly, dcc_biyearly)
# diff_y_b <- monthwise_H_diff(dcc_yearly, dcc_biyearly)

# summary(diff_q_y)
# summary(diff_q_b)
# summary(diff_y_b)

# # Tables and graphs
# print(round(cmp_table, 4))

# comparison_table <- data.frame(
#   Comparison = c("Quarterly vs Yearly", "Quarterly vs Biyearly", "Yearly vs Biyearly"),
#   Frobenius_Distance = as.numeric(comparison_H)
# )

# print(comparison_table)

# diff_summary <- rbind(
#   q_vs_y = summary(diff_q_y),
#   q_vs_b = summary(diff_q_b),
#   y_vs_b = summary(diff_y_b)
# )

# print(round(diff_summary, 3))

# extract_ab <- function(dcc_obj) {
#   data.frame(
#     date = as.Date(names(dcc_obj)),
#     a = sapply(dcc_obj, function(x) x$a),
#     b = sapply(dcc_obj, function(x) x$b),
#     ab = sapply(dcc_obj, function(x) x$a + x$b)
#   )
# }

# df_q <- extract_ab(dcc_quarterly)
# df_y <- extract_ab(dcc_yearly)
# df_b <- extract_ab(dcc_biyearly)

# png("dcc_persistence.png", width = 800, height = 600)
# plot(df_q$date, df_q$ab, type = "l", lwd = 2,
#      main = "DCC Persistence (a+b)",
#      ylab = "a+b", xlab = "")
# lines(df_y$date, df_y$ab, col = 2, lwd = 2)
# lines(df_b$date, df_b$ab, col = 4, lwd = 2)
# legend("topright",
#        legend = c("Quarterly", "Yearly", "Biyearly"),
#        col = c(1,2,4), lwd = 2)
# dev.off()

# png("cov_diff.png", width = 800, height = 800)
# plot(diff_q_y, type = "l", lwd = 2,
#      main = "Monthly Covariance Differences",
#      ylab = "Frobenius Distance", xlab = "Month")
# lines(diff_q_b, col = 2, lwd = 2)
# lines(diff_y_b, col = 4, lwd = 2)
# legend("topright",
#        legend = c("Q vs Y", "Q vs B", "Y vs B"),
#        col = c(1,2,4), lwd = 2)
# dev.off()

# png("mon_by_mon_diff.png", width = 800, height = 800)
# hist(diff_q_y, breaks = 20, col = "gray",
#      main = "Distribution of Covariance Differences (Q vs Y)",
#      xlab = "Frobenius Distance")
# dev.off()

# Adjacency matrix construction 
# @param sigma_hat is the forecasted covariance matrix 
# @return List with partial correlations and the adjacency matrix
adjacency_matrix <- function(sigma_hat, tau = 0.05, ridge = 1e-6,
                             use_correlation = TRUE) {
  sigma_hat <- as.matrix(sigma_hat)
  sigma_hat <- (sigma_hat + t(sigma_hat)) / 2

  if (use_correlation) {
    sigma_hat <- cov_to_cor(sigma_hat, eps = ridge)
    sigma_hat <- make_psd(sigma_hat, eps = ridge)
    diag(sigma_hat) <- 1
  } else {
    sigma_hat <- make_psd(sigma_hat, eps = ridge)
  }

  theta <- solve(sigma_hat + diag(ridge, nrow(sigma_hat)))
  N <- nrow(theta)

  d <- 1 / sqrt(diag(theta))

  partial_corr <- -theta * (d %o% d)
  diag(partial_corr) <- 0

  adjacency_pos <- ifelse(partial_corr > tau, partial_corr, 0)
  adjacency_neg <- ifelse(partial_corr < -tau, -partial_corr, 0)

  return(list(
    partial_corr = partial_corr,
    adjacency_pos = adjacency_pos,
    adjacency_neg = adjacency_neg
  ))
}

# Network Centrality function using EC 
# @param adjacency_pos for positive adjacency matrix
# @param adjacency_neg for negative adjacency matrix 
# @return List with EC+ and EC-
network_centrality <- function(adjacency_pos, adjacency_neg) {
  # positive network
  eig_pos <- eigen(adjacency_pos)
  idx_pos <- which.max(Re(eig_pos$values))
  eigenvector_pos <- Re(eig_pos$vectors[, idx_pos])

  if (sum(eigenvector_pos) < 0) eigenvector_pos <- -eigenvector_pos
  ec_pos <- eigenvector_pos / sum(eigenvector_pos)

  # negative network
  eig_neg <- eigen(adjacency_neg)
  idx_neg <- which.max(Re(eig_neg$values))
  eigenvector_neg <- Re(eig_neg$vectors[, idx_neg])

  if (sum(eigenvector_neg) < 0) eigenvector_neg <- -eigenvector_neg
  ec_neg <- eigenvector_neg / sum(eigenvector_neg)

  return(list(
    ec_pos = ec_pos,
    ec_neg = ec_neg
  ))
}

# Spillover measures to get the systemic risk exposure
# @param adjacency_pos Positive adjacency matrix
# @param adjacency_neg Negative adjacency matrix 
# @param ec_pos EC+
# @param ec_neg EC-
# @param sigma_vec Vector of volatilities
# @return List with spillover_pos and spillover_neg
spillovers <- function(adjacency_pos, adjacency_neg, ec_pos, ec_neg, sigma_vec) {
  # column vector 
  sigma_vec <- as.numeric(sigma_vec)
  # temporary matrix A with zeros on the diagonal for i not j
  temp_adjacency_pos <- adjacency_pos 
  diag(temp_adjacency_pos) <- 0
  temp_adjacency_neg <- adjacency_neg
  diag(temp_adjacency_neg) <- 0

  # positive spillovers 
  init_spillover_pos <- as.numeric(temp_adjacency_pos %*% sigma_vec)
  spillover_pos <- ec_pos * init_spillover_pos

  # negative spillovers
  init_spillover_neg <- as.numeric(temp_adjacency_neg %*% sigma_vec)
  spillover_neg <- ec_neg * init_spillover_neg
  
  # return the list
  return(list(
    spillover_pos = spillover_pos,
    spillover_neg = spillover_neg
  ))
}

# Penalty function g_i(S_i,m) = max(eps, 1 + lambda_pos*S_pos - lambda_neg*S_neg)
network_penalty <- function(spillover_pos, spillover_neg,
                            lambda_pos = 0.1, lambda_neg = 0.1,
                            eps = 1e-4) {
  g <- 1 + lambda_pos * spillover_pos - lambda_neg * spillover_neg
  pmax(eps, g)
}

# Raw network-adjusted weights: \tilde w_i,m = 1 / (sigma_i,m * g_i)
network_adjusted_weights <- function(sigma_vec, spillover_pos, spillover_neg,
                                     lambda_pos = 0.1, lambda_neg = 0.1,
                                     eps = 1e-4) {
  sigma_vec <- as.numeric(sigma_vec)

  g <- network_penalty(
    spillover_pos = spillover_pos,
    spillover_neg = spillover_neg,
    lambda_pos = lambda_pos,
    lambda_neg = lambda_neg,
    eps = eps
  )

  w_tilde <- 1 / (sigma_vec * g)

  list(
    penalty = g,
    w_tilde = w_tilde
  )
}

# Unscaled portfolio return: \tilde f^{NET}_{p,m+1} = sum_i \tilde w_i,m * f_i,m+1
network_portfolio_return <- function(w_tilde, returns_next_month) {
  sum(w_tilde * as.numeric(returns_next_month), na.rm = TRUE)
}

# Scaling constant c so that sd(c * f_net_unscaled) = sd(f_benchmark)
scaling_constant <- function(f_net_unscaled, f_benchmark) {
  sd(f_benchmark, na.rm = TRUE) / sd(f_net_unscaled, na.rm = TRUE)
}

# Final scaled network-managed return
scaled_network_return <- function(f_net_unscaled, c) {
  c * f_net_unscaled
}

run_network_outputs <- function(dcc_list,
                                tau = 0.05,
                                lambda_pos = 0.1,
                                lambda_neg = 0.1,
                                ridge = 1e-6,
                                eps = 1e-4,
                                use_correlation = TRUE,
                                asset_names = NULL) {
  stopifnot(is.list(dcc_list), length(dcc_list) > 0)

  n_periods <- length(dcc_list)
  per_period <- vector("list", n_periods)
  dates <- rep(as.Date(NA), n_periods)

  for (m in seq_len(n_periods)) {
    x <- dcc_list[[m]]
    H_month <- x$H_month
    period_date <- x$refit_date

    if (is.null(H_month)) {
      per_period[[m]] <- list(
        date = period_date,
        H_month = NULL,
        sigma_vec = NULL,
        partial_corr = NULL,
        adjacency_pos = NULL,
        adjacency_neg = NULL,
        ec_pos = NULL,
        ec_neg = NULL,
        spillover_pos = NULL,
        spillover_neg = NULL,
        penalty = NULL,
        w_tilde = NULL
      )
      next
    }

    H_month <- as.matrix(H_month)
    sigma_vec <- sqrt(pmax(diag(H_month), 0))

    current_names <- colnames(H_month)
    if (is.null(current_names)) current_names <- rownames(H_month)
    if (is.null(current_names)) {
      if (!is.null(asset_names)) {
        current_names <- asset_names
      } else {
        current_names <- paste0("Asset_", seq_len(nrow(H_month)))
      }
    }

    adj <- adjacency_matrix(
      sigma_hat = H_month,
      tau = tau,
      ridge = ridge,
      use_correlation = use_correlation
    )

    cent <- network_centrality(
      adjacency_pos = adj$adjacency_pos,
      adjacency_neg = adj$adjacency_neg
    )

    spl <- spillovers(
      adjacency_pos = adj$adjacency_pos,
      adjacency_neg = adj$adjacency_neg,
      ec_pos = cent$ec_pos,
      ec_neg = cent$ec_neg,
      sigma_vec = sigma_vec
    )

    nw <- network_adjusted_weights(
      sigma_vec = sigma_vec,
      spillover_pos = spl$spillover_pos,
      spillover_neg = spl$spillover_neg,
      lambda_pos = lambda_pos,
      lambda_neg = lambda_neg,
      eps = eps
    )

    per_period[[m]] <- list(
      date = period_date,
      H_month = H_month,
      sigma_vec = setNames(as.numeric(sigma_vec), current_names),

      partial_corr = `dimnames<-`(adj$partial_corr, list(current_names, current_names)),
      adjacency_pos = `dimnames<-`(adj$adjacency_pos, list(current_names, current_names)),
      adjacency_neg = `dimnames<-`(adj$adjacency_neg, list(current_names, current_names)),

      ec_pos = setNames(as.numeric(cent$ec_pos), current_names),
      ec_neg = setNames(as.numeric(cent$ec_neg), current_names),

      spillover_pos = setNames(as.numeric(spl$spillover_pos), current_names),
      spillover_neg = setNames(as.numeric(spl$spillover_neg), current_names),

      penalty = setNames(as.numeric(nw$penalty), current_names),
      w_tilde = setNames(as.numeric(nw$w_tilde), current_names)
    )

    dates[m] <- as.Date(period_date)
  }

  ok <- sapply(per_period, function(z) !is.null(z$w_tilde))
  if (!any(ok)) stop("No valid network outputs produced.")

  w_tilde_mat <- do.call(rbind, lapply(per_period[ok], function(z) unname(z$w_tilde)))
  penalty_mat <- do.call(rbind, lapply(per_period[ok], function(z) unname(z$penalty)))
  spillover_pos_mat <- do.call(rbind, lapply(per_period[ok], function(z) unname(z$spillover_pos)))
  spillover_neg_mat <- do.call(rbind, lapply(per_period[ok], function(z) unname(z$spillover_neg)))

  first_ok <- which(ok)[1]
  final_names <- names(per_period[[first_ok]]$w_tilde)

  colnames(w_tilde_mat) <- final_names
  colnames(penalty_mat) <- final_names
  colnames(spillover_pos_mat) <- final_names
  colnames(spillover_neg_mat) <- final_names

  row_ids <- as.character(dates[ok])
  rownames(w_tilde_mat) <- row_ids
  rownames(penalty_mat) <- row_ids
  rownames(spillover_pos_mat) <- row_ids
  rownames(spillover_neg_mat) <- row_ids

  summary_df <- data.frame(
    period = which(ok),
    date = dates[ok],
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_df,
    w_tilde = w_tilde_mat,
    penalty = penalty_mat,
    spillover_pos = spillover_pos_mat,
    spillover_neg = spillover_neg_mat,
    per_period = per_period
  )
}

extract_realized_forecast_returns <- function(dcc_list, returns_df) {
  returns_dates <- as.Date(returns_df$date)
  returns_mat <- as.matrix(returns_df[, -1, drop = FALSE])

  out_list <- vector("list", length(dcc_list))
  out_dates <- rep(as.Date(NA), length(dcc_list))

  for (m in seq_along(dcc_list)) {
    x <- dcc_list[[m]]

    if (is.null(x$forecast_dates) || length(x$forecast_dates) == 0) {
      out_list[[m]] <- NULL
      next
    }

    fc_dates <- as.Date(x$forecast_dates)
    idx <- match(fc_dates, returns_dates)

    if (any(is.na(idx))) {
      stop(sprintf("Could not match forecast dates to returns_df in period %d.", m))
    }

    ret_block <- returns_mat[idx, , drop = FALSE]

    realized_vec <- apply(ret_block, 2, function(z) prod(1 + z, na.rm = TRUE) - 1)

    out_list[[m]] <- realized_vec
    out_dates[m] <- as.Date(x$refit_date)
  }

  ok <- !sapply(out_list, is.null)
  realized_mat <- do.call(rbind, out_list[ok])

  colnames(realized_mat) <- colnames(returns_df[, -1, drop = FALSE])
  rownames(realized_mat) <- as.character(out_dates[ok])

  list(
    summary = data.frame(
      period = which(ok),
      date = out_dates[ok],
      stringsAsFactors = FALSE
    ),
    realized_returns = realized_mat,
    per_period = out_list
  )
}

compute_network_returns <- function(net_obj, realized_obj) {
  common_dates <- intersect(rownames(net_obj$w_tilde), rownames(realized_obj$realized_returns))

  W <- net_obj$w_tilde[common_dates, , drop = FALSE]
  R <- realized_obj$realized_returns[common_dates, , drop = FALSE]

  if (!identical(colnames(W), colnames(R))) {
    stop("Column names of weights and realized returns do not match.")
  }

  data.frame(
    date = as.Date(common_dates),
    network_return_unscaled = rowSums(W * R, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

compute_equal_weight_benchmark <- function(realized_obj, common_dates = NULL) {
  R <- realized_obj$realized_returns

  if (!is.null(common_dates)) {
    R <- R[common_dates, , drop = FALSE]
  }

  ew_ret <- rowMeans(R, na.rm = TRUE)

  data.frame(
    date = as.Date(rownames(R)),
    EW_benchmark = as.numeric(ew_ret),
    stringsAsFactors = FALSE
  )
}

# # =================================================
# # Run for current test case
# # first 1100 obs, first 8 factors
# # =================================================
# factor_names_8 <- colnames(test_data_cmp[, -1, drop = FALSE])

# net_biyearly <- run_network_outputs(
#   dcc_list = dcc_biyearly,
#   # tau = sqrt(0.05),
#   tau = 0.05, 
#   lambda_pos = 0.1,
#   lambda_neg = 0.1,
#   asset_names = factor_names_8
# )

# realized_biyearly <- extract_realized_forecast_returns(
#   dcc_list = dcc_biyearly,
#   returns_df = test_data_cmp
# )

# # force network output column names to match realized returns exactly
# colnames(net_biyearly$w_tilde) <- colnames(realized_biyearly$realized_returns)
# colnames(net_biyearly$penalty) <- colnames(realized_biyearly$realized_returns)
# colnames(net_biyearly$spillover_pos) <- colnames(realized_biyearly$realized_returns)
# colnames(net_biyearly$spillover_neg) <- colnames(realized_biyearly$realized_returns)

# network_returns_biyearly <- compute_network_returns(
#   net_obj = net_biyearly,
#   realized_obj = realized_biyearly
# )

# benchmark_biyearly <- compute_equal_weight_benchmark(
#   realized_obj = realized_biyearly,
#   common_dates = as.character(network_returns_biyearly$date)
# )

# # optional combined object for later performance code
# network_vs_benchmark_biyearly <- network_returns_biyearly |>
#   left_join(benchmark_biyearly, by = "date")

# # inspect
# net_biyearly$summary
# net_biyearly$w_tilde
# net_biyearly$penalty
# net_biyearly$spillover_pos
# net_biyearly$spillover_neg

# realized_biyearly$realized_returns

# network_returns_biyearly
# benchmark_biyearly
# network_vs_benchmark_biyearly

# =================================================
# Run network pipeline for current full-sample setup
# =================================================

factor_names_all <- colnames(managed_portfolios[, -1, drop = FALSE])

net_all <- run_network_outputs(
  dcc_list = dcc_all_full,
  tau = 0.05,
  lambda_pos = 0.1,
  lambda_neg = 0.1,
  asset_names = factor_names_all
)

realized_all <- extract_realized_forecast_returns(
  dcc_list = dcc_all_full,
  returns_df = managed_portfolios
)

# force network output column names to match realized returns exactly
colnames(net_all$w_tilde)       <- colnames(realized_all$realized_returns)
colnames(net_all$penalty)       <- colnames(realized_all$realized_returns)
colnames(net_all$spillover_pos) <- colnames(realized_all$realized_returns)
colnames(net_all$spillover_neg) <- colnames(realized_all$realized_returns)

network_returns_all <- compute_network_returns(
  net_obj = net_all,
  realized_obj = realized_all
)

benchmark_all <- compute_equal_weight_benchmark(
  realized_obj = realized_all,
  common_dates = as.character(network_returns_all$date)
)

network_vs_benchmark_all <- network_returns_all |>
  left_join(benchmark_all, by = "date")

# inspect
# net_all$summary
# net_all$w_tilde
# net_all$penalty
# net_all$spillover_pos
# net_all$spillover_neg

# realized_all$realized_returns

# network_returns_all
# benchmark_all
# network_vs_benchmark_all

# Factor labels (nicer names in plots)
factor_labels <- c(
  mkt_excess = "MKT",
  smb = "SMB",
  hml = "HML",
  rmw = "RMW",
  cma = "CMA",
  mom = "MOM",
  agric = "AGRIC",
  food = "FOOD"
)

# Build graph
build_factor_graph <- function(adj_matrix) {
  A <- as.matrix(adj_matrix)
  diag(A) <- 0

  graph_from_adjacency_matrix(
    A,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
}

# Plot single network
plot_factor_network <- function(net_obj, date, type = "positive") {

  idx <- which(as.character(net_obj$summary$date) == as.character(as.Date(date)))
  period_id <- net_obj$summary$period[idx]
  pp <- net_obj$per_period[[period_id]]

  if (type == "positive") {
    A <- pp$adjacency_pos
    edge_col <- "darkgreen"
  } else {
    A <- pp$adjacency_neg
    edge_col <- "red"
  }

  g <- build_factor_graph(A)

  # labels
  V(g)$label <- factor_labels[V(g)$name]

  # node size (centrality)
  cent <- if (type == "positive") pp$ec_pos else pp$ec_neg
  cent <- cent[V(g)$name]

  V(g)$size <- 10 + 30 * (cent / max(cent))

  # edge width
  E(g)$width <- 1 + 6 * (E(g)$weight / max(E(g)$weight))

  V(g)$color <- "lightblue"
  E(g)$color <- edge_col

  plot(
    g,
    layout = layout_in_circle(g),
    main = paste(type, "network -", date)
  )
}

# Plot positive + negative side by side
plot_network_pair <- function(net_obj, date) {
  par(mfrow = c(1, 2))

  plot_factor_network(net_obj, date, "positive")
  plot_factor_network(net_obj, date, "negative")

  par(mfrow = c(1, 1))
}

# Plot sequence
plot_network_sequence <- function(net_obj, dates, type = "positive") {
  n <- length(dates)
  par(mfrow = c(ceiling(n / 2), 2))

  for (d in dates) {
    plot_factor_network(net_obj, d, type)
  }

  par(mfrow = c(1, 1))
}

# Heatmap of partial correlations
plot_heatmap <- function(net_obj, date) {

  idx <- which(as.character(net_obj$summary$date) == as.character(as.Date(date)))
  period_id <- net_obj$summary$period[idx]
  P <- net_obj$per_period[[period_id]]$partial_corr

  lab <- factor_labels[colnames(P)]
  colnames(P) <- lab
  rownames(P) <- lab

  image(
    1:nrow(P), 1:ncol(P), t(P[nrow(P):1, ]),
    col = colorRampPalette(c("red", "white", "green"))(100),
    axes = FALSE,
    main = paste("Partial correlations -", date)
  )

  axis(1, at = 1:ncol(P), labels = colnames(P), las = 2)
  axis(2, at = 1:nrow(P), labels = rev(rownames(P)), las = 2)
}

# 1) Positive + Negative network (side-by-side)
png("network_pair_1974_09_30.png", width = 1200, height = 600)
plot_network_pair(net_biyearly, "1974-09-30")
dev.off()

# 2) Positive network only
png("network_positive_1974_09_30.png", width = 800, height = 600)
plot_factor_network(net_biyearly, "1974-09-30", "positive")
dev.off()

# 3) Negative network only
png("network_negative_1974_09_30.png", width = 800, height = 600)
plot_factor_network(net_biyearly, "1974-09-30", "negative")
dev.off()

# 4) Sequence of first 6 months (positive networks)
png("network_sequence_positive.png", width = 1400, height = 900)
plot_network_sequence(
  net_biyearly,
  net_biyearly$summary$date[1:6],
  "positive"
)
dev.off()

# 5) Heatmap
png("heatmap_1974_09_30.png", width = 800, height = 700)
plot_heatmap(net_biyearly, "1974-09-30")
dev.off()

# =================================================
# Benchmarks
# =================================================
# Setting here to test the script for less managed portfolios
#managed_portfolios <- managed_portfolios[,1:7]

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

#Scaling b for target volatility
vol_mve_in_sample <- sd(as.matrix(managed_portfolios_est[,-1]) %*% b) *sqrt(252)
b <- b * 0.1 /vol_mve_in_sample

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
    rv = sum((MVE_return-mean(MVE_return))^2),
    .groups = "drop"
  ) |> mutate(
    rv_lag = lag(rv)
  ) |> select(-rv)

MVE_returns_df <- MVE_returns_df |> group_by(month) |> 
  summarise(
  monthly_return = prod(1 + MVE_return) - 1
)

# Scaling the portfolio by its lagged realized variance
MVE_returns_df <- MVE_returns_df |> left_join(rv_lag_MVE, by = "month")
MVE_returns_df <- MVE_returns_df |> 
  mutate(
    MVE_scaled = monthly_return / rv_lag
  )
MVE_returns_df <- MVE_returns_df |> select(-c("rv_lag", "n_of_days"))

# Take out the first month, because we do not have the realized variance for the period 0
MVE_returns_df <- MVE_returns_df |> filter(month >= start_date_estimation + months(1))

# Scaling by c
c = sd(MVE_returns_df$monthly_return)/sd(MVE_returns_df$MVE_scaled)
MVE_returns_df$MVE_strategy_return <- MVE_returns_df$MVE_scaled*c
MVE_returns_df <- MVE_returns_df |> select(-c(monthly_return, MVE_scaled))


### Equally weighted buy-and-hold 
no_of_factors <- dim(managed_portfolios[,-1])[2]
b_EW <- as.matrix(rep(1/no_of_factors, no_of_factors))

EW_returns <- as.matrix(managed_portfolios[,-1]) %*% b_EW
EW_returns_df <- tibble(
  date = managed_portfolios$date, 
  month = floor_date(date, "month"),
  EW_return = as.numeric(EW_returns)
) 

EW_returns_df <- EW_returns_df |> group_by(month) |>
  summarise(
    monthly_return = prod(1 + EW_return) - 1
  ) |> filter(month >= start_date_estimation + months(1))


# Changing data frames only for out of sample periods
MVE_returns_df <- MVE_returns_df |> filter(month > end_date_estimation)
EW_returns_df <- EW_returns_df |> filter(month > end_date_estimation)

## This is just for verifying the performance of the benchmarks
# Calculating the Sharpe Ratios
sharpe_ratios_benchmarks <- tibble(
  Benchmark = c("EW Buy-and-Hold", "MVE"),
  SR = c(round(compute_SR(EW_returns_df$monthly_return),4)*sqrt(12), 
         round(compute_SR(MVE_returns_df$MVE_strategy_return),4)*sqrt(12))
)


# Creating data frame for letter regressions
benchmarks_returns <- EW_returns_df |> left_join(MVE_returns_df, by = "month") 
colnames(benchmarks_returns) <- c("month", "EW", "MVE")

# Functions for getting the performance measures for network analysis 
# calculate scaling constant c to match volatility of BH 
vol_target <- sd(network_vs_benchmark_all$EW_benchmark, na.rm = TRUE)
vol_strat_unscaled <- sd(network_vs_benchmark_all$network_return_unscaled, na.rm = TRUE)
c_scaling <- vol_target / vol_strat_unscaled

network_vs_benchmark_all <- network_vs_benchmark_all %>%
  mutate(net_strategy_return = network_return_unscaled * c_scaling)

# ====================================
# Generating figures 
# ====================================
# Cumulative wealth dataframe
# Prepare the data by joining NET returns to the benchmarks (UNCOMMENT THIS WHEN RUNNING FULL SAMPLE)
# Build monthly network returns with aligned dates
network_monthly <- network_vs_benchmark_all %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    net_strategy_return = first(net_strategy_return),
    .groups = "drop"
  )

# Join to benchmark monthly returns
combined_returns <- benchmarks_returns %>%
  left_join(network_monthly, by = "month") %>%
  arrange(month)

# Build cumulative wealth safely
first_net_idx <- which(!is.na(combined_returns$net_strategy_return))[1]

cumulative_wealth_df <- combined_returns %>%
  mutate(
    EW = cumprod(1 + EW),
    MVE = cumprod(1 + MVE),
    NET = NA_real_
  )

if (!is.na(first_net_idx)) {
  cumulative_wealth_df$NET[first_net_idx:nrow(cumulative_wealth_df)] <-
    cumprod(1 + cumulative_wealth_df$net_strategy_return[first_net_idx:nrow(cumulative_wealth_df)])
}

png("cumulative_wealth.png", width = 900, height = 600)

axis <- par(lab = c(20, 8, 5))
max_y <- max(cumulative_wealth_df$NET, cumulative_wealth_df$MVE, na.rm = TRUE)

plot(x = cumulative_wealth_df$month, y = cumulative_wealth_df$EW,
     xlab = "Date", ylab = "Cumulative return",
     type = "l", col = "black", lwd = 1, lty = 1, ylim = c(0, max_y))

y_ticks <- pretty(c(0, max_y))
abline(h = y_ticks, col = "grey85", lty = 1)

lines(x = cumulative_wealth_df$month, y = cumulative_wealth_df$EW, col = "black", lwd = 1)
lines(x = cumulative_wealth_df$month, y = cumulative_wealth_df$MVE, col = "blue", lty = 2, lwd = 1.5)
lines(x = cumulative_wealth_df$month, y = cumulative_wealth_df$NET, col = "red", lty = 3, lwd = 2)

legend("topleft", legend = c("BH", "MVE", "NET"),
       col = c("black", "blue", "red"), lty = c(1, 2, 3),
       lwd = 2, bty = "n", cex = 0.8)

dev.off()

# Rolling Sharpe Ratio
rolling_SR_df <- data.frame(
  date = combined_returns$month[12:nrow(combined_returns)],
  SR_EW = NA_real_,
  SR_MVE = NA_real_,
  SR_NET = NA_real_
)

for(i in 12:nrow(combined_returns)){
  rolling_SR_df$SR_EW[i-11]  <- mean(combined_returns$EW[(i-11):i]) / sd(combined_returns$EW[(i-11):i]) * sqrt(12)
  rolling_SR_df$SR_MVE[i-11] <- mean(combined_returns$MVE[(i-11):i]) / sd(combined_returns$MVE[(i-11):i]) * sqrt(12)
  rolling_SR_df$SR_NET[i-11] <- mean(combined_returns$net_strategy_return[(i-11):i]) / sd(combined_returns$net_strategy_return[(i-11):i]) * sqrt(12) # Added this
}

png("rolling_sharpe.png", width = 900, height = 600)

plot(x = rolling_SR_df$date, y = rolling_SR_df$SR_EW,
     type = "l", xlab = "Date", ylab = "Rolling SR",
     lty = 1, col = "black",
     ylim = range(c(rolling_SR_df$SR_EW,
                    rolling_SR_df$SR_MVE,
                    rolling_SR_df$SR_NET), na.rm = TRUE))

y_ticks <- pretty(rolling_SR_df$SR_EW)
abline(h = y_ticks, col = "grey85", lty = 1)

lines(x = rolling_SR_df$date, y = rolling_SR_df$SR_EW, col = "black", lty = 1)
lines(x = rolling_SR_df$date, y = rolling_SR_df$SR_MVE, col = "blue", lty = 2)
lines(x = rolling_SR_df$date, y = rolling_SR_df$SR_NET, col = "red", lty = 3, lwd = 2)

legend("topright", legend = c("BH", "MVE", "NET"),
       col = c("black", "blue", "red"),
       lty = c(1, 2, 3), bty = "n", lwd = 2, cex = 0.8)

dev.off()
# Drawdown figure
drawdown_df <- cumulative_wealth_df |>
  arrange(month) |>
  mutate(
    EW_drawdown  = EW / cummax(EW) - 1,
    MVE_drawdown = MVE / cummax(MVE) - 1,
    NET_drawdown = NA_real_
  )

net_valid_idx <- which(!is.na(drawdown_df$NET))

if (length(net_valid_idx) > 0) {
  first_net <- min(net_valid_idx)
  net_path <- drawdown_df$NET[first_net:nrow(drawdown_df)]
  drawdown_df$NET_drawdown[first_net:nrow(drawdown_df)] <- net_path / cummax(net_path) - 1
}

png("drawdown.png", width = 900, height = 600)

plot(x = drawdown_df$month, y = drawdown_df$EW_drawdown,
     type = "l", ylim = c(-0.6, 0),
     xlab = "Date", ylab = "Drawdown")

y_ticks <- pretty(drawdown_df$EW_drawdown)
abline(h = y_ticks, col = "grey85", lty = 1)

lines(x = drawdown_df$month, y = drawdown_df$EW_drawdown, col = "black", lty = 1)
lines(x = drawdown_df$month, y = drawdown_df$MVE_drawdown, col = "blue", lty = 2)
lines(x = drawdown_df$month, y = drawdown_df$NET_drawdown, col = "red", lty = 3, lwd = 2)

legend("bottomright", legend = c("BH", "MVE", "NET"),
       col = c("black", "blue", "red"),
       lty = c(1, 2, 3), bty = "n", lwd = 2, cex = 0.8)

dev.off()

sum(is.na(cumulative_wealth_df$NET))
sum(is.na(cumulative_wealth_df$MVE))

sum(is.infinite(cumulative_wealth_df$NET))
sum(is.infinite(cumulative_wealth_df$MVE))

which(is.na(cumulative_wealth_df$NET))[1:10]

cumulative_wealth_df[which(is.na(cumulative_wealth_df$NET)), ]
# ================================================
# Performance evaluation and alpha testing
# ================================================
#scale = 12 for monthly data -> annualising the metrics
# scale = 252 for daily 
annualized_mean <- function(r, scale = 12) {
  mean(r, na.rm = TRUE) * scale
}

annualized_vol <- function(r, scale = 12) {
  sd(r, na.rm = TRUE) * sqrt(scale)
}

compute_SR <- function(r, scale = 252) {
  mean(r, na.rm = TRUE) / sd(r, na.rm = TRUE) * sqrt(scale)
}

max_drawdown <- function(r) {
  wealth <- cumprod(1 + r)
  drawdown <- wealth / cummax(wealth) - 1
  min(drawdown, na.rm = TRUE)
}

compute_turnover <- function(weights_df) {
  #this is assuming the first column is a date or something 
  W <- as.matrix(weights_df[,-1])
  turnover <- c(NA, rowSums(abs(W[-1,] - W[-nrow(W),]), na.rm = TRUE))
  return(turnover)
}

compute_net_returns <- function(r, turnover, cost = 0.001) {
  r - cost * turnover
}

performance_summary <- function(r, turnover = NULL, scale = 12) {
  data.frame(
    Mean = round(annualized_mean(r, scale),4) * 100,
    Volatility = round(annualized_vol(r, scale),4) * 100,
    Sharpe = round(compute_SR(r, scale),2),
    MaxDrawdown = round(max_drawdown(r),4) * 100,
    AvgTurnover = ifelse(is.null(turnover), NA, round(mean(turnover, na.rm = TRUE),4)),
    NetMean = ifelse(is.null(turnover), NA, round(annualized_mean(compute_net_returns(r, turnover),4), scale)),
    NetSharpe = ifelse(is.null(turnover), NA, round(compute_SR(compute_net_returns(r, turnover),4), scale))
  )
}
perf_EW  <- performance_summary(benchmarks_returns$EW)
perf_MVE <- performance_summary(benchmarks_returns$MVE)
perf_NET <- performance_summary(network_vs_benchmark_biyearly$net_strategy_return)

#just to make a nice table with the different strategy
performance_table <- bind_rows(
  `EW Buy-and-Hold` = perf_EW,
  `MVE Vol-Managed` = perf_MVE,
  `NET Strategy`    = perf_NET,
  .id = "Strategy"
)

print(performance_table)

alpha_test <- function(strategy_ret, benchmark_ret, lag = 3) {
  df <- data.frame(
    y = strategy_ret,
    x = benchmark_ret
  ) %>%
    drop_na()
  
  fit <- lm(y ~ x, data = df)
  nw_se <- NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = TRUE)
  test <- coeftest(fit, vcov. = nw_se)
  
  data.frame(
    alpha = test["(Intercept)", "Estimate"]* 12,
    alpha_se = test["(Intercept)", "Std. Error"]* 12 ,
    alpha_t = test["(Intercept)", "t value"] ,
    alpha_p = test["(Intercept)", "Pr(>|t|)"] ,
    beta = test["x", "Estimate"] * 12
  )
}
alpha_vs_EW <- alpha_test(
  strategy_ret = benchmarks_returns$MVE,
  benchmark_ret = benchmarks_returns$EW
)
print(alpha_vs_EW)

alpha_NET_vs_BH <- alpha_test(
  strategy_ret = network_vs_benchmark_biyearly$net_strategy_return,
  benchmark_ret = network_vs_benchmark_biyearly$EW_benchmark
)

print(alpha_NET_vs_BH)

# Make sure MVE is in the comparison table
network_vs_benchmark_biyearly <- network_vs_benchmark_biyearly %>%
  left_join(benchmarks_returns %>% select(month, MVE), by = c("date" = "month"))

alpha_NET_vs_MVE <- alpha_test(
  strategy_ret = network_vs_benchmark_biyearly$net_strategy_return,
  benchmark_ret = network_vs_benchmark_biyearly$MVE 
)

print(alpha_NET_vs_MVE)

sharpe_difference <- function(r1, r2, scale = 12) {
  compute_SR(r1, scale) - compute_SR(r2, scale)
}
sr_diff <- sharpe_difference(
  benchmarks_returns$MVE,
  benchmarks_returns$EW
)
print(sr_diff)
# Ledoit-Wolf Sharpe ratio difference test
lw_sharpe_test <- sharpeTesting(
  x = benchmarks_returns$MVE_strategy_return,
  y = benchmarks_returns$EW_return,
  control = list(
    type = 2,
    hac = TRUE,
    nBoot = 999,
    bBoot = 1
  )
)

lw_summary <- tibble(
  Strategy_1 = "MVE",
  Strategy_2 = "EW Buy-and-Hold",
  SR_1 = lw_sharpe_test$sharpe[1],
  SR_2 = lw_sharpe_test$sharpe[2],
  SR_Difference = lw_sharpe_test$dsharpe,
  t_stat = lw_sharpe_test$tstat,
  p_value = lw_sharpe_test$pval
)

print(lw_summary)

# Difference between NET and Buy-and-Hold 
sr_diff_net_ew <- sharpe_difference(
  network_vs_benchmark_biyearly$net_strategy_return,
  network_vs_benchmark_biyearly$EW_benchmark
)
print(sr_diff_net_ew)

# Difference between NET and MVE 
sr_diff_net_mve <- sharpe_difference(
  network_vs_benchmark_biyearly$net_strategy_return,
  network_vs_benchmark_biyearly$MVE
)
print(sr_diff_net_mve)

# ========================================
# Table creation
# ========================================





#' Function creates the Panel A of the main results table, 
#' it contains gross mean, volatility, SR, and alphas for 
#' both benchmarks
#'
#'@param returns_df Data frame with date, and returns corresponding to:
#' equally-weighted-buy and hold (EW), mean-variance efficient (MVE) and
#' netowrk strategy (NET)
#'
#'@return Panel A of the main results table
#'
main_results_table_panel_A <- function(returns_df){
  
  if(class(returns_df[[1]]) == "Date"){
    returns_df <- returns_df[,-1]
  }
  
  # Computing the input for the table 
  gross_mean_return <- colMeans(returns_df) * 12 * 100
  gross_volatility <- sapply(returns_df, sd) * sqrt(12) * 100
  SR <- gross_mean_return/gross_volatility
  alpha_BH <- c(0,0,0)
  alpha_MVE <- c(0,0,0)
  
  results <- data.frame(
    Strategy = c("BH", "MVE", "NET"),
    Gross_mean_return = gross_mean_return,
    Gross_volatility = gross_volatility,
    Sharpe_ratio = SR,
    Alpha_BH <- alpha_BH,
    Alpha_MVE <- alpha_MVE
  )
  
  
  return(results)
}


#'
#'Function creates the sub-period analysis table,
#'it contians mean gross return, gross volatility, SR, 
#'alphas and t-statistics
#'
#'@param returns_df Data frame with date, and returns corresponding to:
#' equally-weighted-buy and hold (EW), mean-variance efficient (MVE) and
#' netowrk strategy (NET)
#'
#'@return Table
#'
sub_period_analysis_table <- function(returns_df){
  
  # Setting dates for first period
  start_date_period_1 <- as.Date("1973-01-01")
  end_date_period_1 <- as.Date("2009-12-13")
  returns_df_period_1 <- returns_df |> filter (date >= start_date_period_1 & date <= end_date_period_1)
  
  # Setting dates for second period
  start_date_period_2 <- as.Date("2010-01-01")
  end_date_period_2 <- as.Date("2017-12-31")
  returns_df_period_2 <- returns_df |> filter (date >= start_date_period_2 & date <= end_date_period_2)
  
  # Setting dates for third period
  start_date_period_3 <- as.Date("2018-01-01")
  end_date_period_3 <- as.Date("2026-01-31")
  returns_df_period_3 <- returns_df |> filter (date >= start_date_period_3 & date <= end_date_period_3)
  
  if(class(returns_df[[1]]) == "Date"){
    returns_df <- returns_df[,-1]
  }
  
  # Calculating gross mean returns
  gross_mean_return_1 <- colMeans(returns_df_period_1) * 12 * 100
  gross_mean_return_2 <- colMeans(returns_df_period_2) * 12 * 100
  gross_mean_return_3 <- colMeans(returns_df_period_3) * 12 * 100
  
  # Calculating gross volaitlities
  gross_volatility_1 <- sapply(returns_df_period_1, sd) * sqrt(12) * 100
  gross_volatility_2 <- sapply(returns_df_period_2, sd) * sqrt(12) * 100
  gross_volatility_3 <- sapply(returns_df_period_3, sd) * sqrt(12) * 100
  
  # Calculating SR's
  SR_1 <- gross_mean_return_1/gross_volatility_1
  SR_2 <- gross_mean_return_2/gross_volatility_2
  SR_3 <- gross_mean_return_3/gross_volatility_3
  
  # Calcualting alphas BH
  alpha_BH_1 <- c(0,0,0)
  alpha_BH_2 <- c(0,0,0)
  alpha_BH_3 <- c(0,0,0)
  
  # Calcualting alphas MVE
  alpha_MVE_1 <- c(0,0,0)
  alpha_MVE_2 <- c(0,0,0)
  alpha_MVE_3 <- c(0,0,0)
  
  
  results <- data.frame(
    Period = c("1973-2010", "", "", "2010-2018", "", "", "2018-2026", "", ""),
    Strategy = c("BH", "MVE", "NET", "BH", "MVE", "NET", "BH", "MVE", "NET"),
    Mean = rbind(gross_mean_return_1, gross_mean_return_2, gross_mean_return_3),
    Volatility = rbind(gross_volatility_1, gross_volatility_2, gross_volatility_3),
    Sharpe_ratio = rbind(SR_1, SR_2, SR_3), 
    Alpha_EW = rbind(alpha_BH_1, alpha_BH_2, alpha_BH_3), 
    Alpha_MVE = rbind(alpha_MVE_1, alpha_MVE_2, alpha_MVE_3)
  )
  
}