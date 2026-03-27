##### Uncomment if packages not installed
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
library(igraph)

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
# =================================================

#' Compute the MVE static weights
#' 
#' 
#' @param mu Vector of means
#' @param sigma Variance-covariance matrix of returns
#' @return Solved weights b
compute_MVE_weights <- function(mu, sigma){
  b <- solve(sigma, mu)
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

# Keeping track of current month and the forecasting month
make_month_index <- function(dates, min_obs = 252) {
  dates <- as.Date(dates)
  n <- length(dates)

  if (n == 0) stop("dates is empty.")

  month_id <- format(dates, "%Y-%m")

  # End-of-month row index for each observed month
  month_end_idx <- tapply(seq_len(n), month_id, max)
  month_end_idx <- as.integer(month_end_idx)

  # Creating a list of monthly forecasting blocks:
  # for each month-end t, forecasts the rows belonging to next month
  out <- vector("list", length(month_end_idx) - 1)
  out_pos <- 1

  for (k in seq_len(length(month_end_idx) - 1)) {
    t_end <- month_end_idx[k]
    next_end <- month_end_idx[k + 1]

    insample_idx <- seq_len(t_end)
    forecast_idx <- (t_end + 1):next_end

    if (length(insample_idx) < min_obs) next
    if (length(forecast_idx) == 0) next

    out[[out_pos]] <- list(
      refit_month = names(month_end_idx)[k],
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
  S_target <- make_psd(S_target, eps = eps)
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
    # Enforce symmetry 
    Q_t <- (Q_t + t(Q_t)) / 2
    R_t <- normalize_Q_to_R(Q_t)

    Q_list[[t]] <- Q_t
    R_list[[t]] <- R_t
    Q_prev <- Q_t
  }

  list(Q = Q_list, R = R_list, Q_T = Q_prev)
}

# Computing negative loglikelihood of the DCC model given (a,b)
dcc_negloglik <- function(par, Z, S, penalty = 1e12) {
  # Parameter vector
  a <- par[1]
  b <- par[2]

  # DCC constraints
  if (!is.finite(a) || !is.finite(b) || a < 0 || b < 0 || (a + b) >= 0.999) {
    return(penalty)
  }

  # Just some matrix formatting, more for clarity
  Z <- as.matrix(Z)
  Tn <- nrow(Z)
  N <- ncol(Z)

  # Given (a,b), first run the entire DCC filter to construct all correlation matrices
  filt <- dcc_filter(Z, a = a, b = b, S = S)
  nll <- 0


  for (t in seq_len(Tn)) {
    zt <- Z[t, ]
    if (any(!is.finite(zt))) next

    # getting model implied R_t, and then make it positive semidefinite
    Rt <- filt$R[[t]]
    Rt <- make_psd(Rt)

    detR <- determinant(Rt, logarithm = TRUE)
    if (!is.finite(detR$modulus)) return(penalty)

    invR <- tryCatch(solve(Rt), error = function(e) NULL)
    if (is.null(invR)) return(penalty)

    quad <- drop(t(zt) %*% invR %*% zt)
    if (!is.finite(quad)) return(penalty)

    nll <- nll + as.numeric(detR$modulus) + quad
  }

  # One-half factor as in classic Gaussian Likelihood
  0.5 * nll
}

# Function to estimate DCC parameters (a,b) from in-sample standardised residuals
estimate_dcc <- function(Z_insample, S_target = NULL, 
start_par = c(0.03, 0.95)) { # Need to justify starting values
  Z_insample <- as.matrix(Z_insample)

  # Drop rows with any missing values for DCC estimation
  keep <- complete.cases(Z_insample)
  Z_use <- Z_insample[keep, , drop = FALSE]

  if (nrow(Z_use) < 50) {
    stop("Too few complete observations to estimate DCC.")
  }

  # Choosing long-run target S
  if (is.null(S_target)) {
    S_target <- safe_cor(Z_use)
  } else {
    S_target <- make_psd(S_target)
    S_target <- cov_to_cor(S_target)
  }

  # Optimisation function to minimise log likelihood
  # We use L-BFGS to minmise computation, and its also better for few parameters
  opt <- optim(
    par = start_par,
    fn = dcc_negloglik,
    Z = Z_use,
    S = S_target,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.999)
  )

  # Estimated Parameters
  a_hat <- opt$par[1]
  b_hat <- opt$par[2]

  # Enforce stationarity softly if optimiser lands near the boundary
  if ((a_hat + b_hat) >= 0.999) {
    s <- a_hat + b_hat
    a_hat <- a_hat * 0.999 / s
    b_hat <- b_hat * 0.999 / s
  }

  # Gives in-sample state Q_t, which we need to forecast on
  filt <- dcc_filter(Z_use, a = a_hat, b = b_hat, S = S_target)

  list(
    a = a_hat,
    b = b_hat,
    S = S_target,
    Q_T = filt$Q_T,
    opt = opt,
    n_obs = nrow(Z_use)
  )
}

# Function that uses DCC estimation to forecast multiple steps ahead correlations
forecast_dcc_correlations <- function(Q_T, S, a, b, h) {
  phi <- a + b # Persistence term
  out <- vector("list", h)

  Q_prev <- Q_T
  # forecasting from the last in-sample state and move forward day by day
  for (k in seq_len(h)) {
    # Expected future DCC recursion
    Q_fc <- S + phi * (Q_prev - S)
    Q_fc <- (Q_fc + t(Q_fc)) / 2
    R_fc <- normalize_Q_to_R(Q_fc)

    out[[k]] <- list(Q = Q_fc, R = R_fc)
    Q_prev <- Q_fc
  }

  out
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
    R_t <- R_list[[t]]$R
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

# Master function
# loops over month-end refit dates and does the full monthly DCC forecasting exercise
# estimate DCC on expanding in-sample data
# forecast next-month daily correlations
# combine with next-month daily volatility forecasts
# aggregate into monthly covariance forecasts
run_dcc_monthly <- function(Z_block, SIGMA_block, dates_block,
                            min_corr_window = 252, S_builder = NULL,
                            trace = TRUE) {
  Z_block <- as.matrix(Z_block)
  SIGMA_block <- as.matrix(SIGMA_block)
  dates_block <- as.Date(dates_block)

  if (!identical(dim(Z_block), dim(SIGMA_block))) {
    stop("Z_block and SIGMA_block must have identical dimensions.")
  }

  if (nrow(Z_block) != length(dates_block)) {
    stop("dates_block length must equal nrow(Z_block).")
  }

  # monthly expanding-window forecasting schedule
  month_map <- make_month_index(dates_block, min_obs = min_corr_window)

  if (length(month_map) == 0) {
    stop("No eligible monthly DCC forecasting blocks found.")
  }

  out <- vector("list", length(month_map))

  for (i in seq_along(month_map)) {
    blk <- month_map[[i]]
    
    # Diagnostic measure since DCC fitting can take time
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

    Z_insample <- Z_block[blk$insample_idx, , drop = FALSE] # DCC estimation sample
    SIGMA_future <- SIGMA_block[blk$forecast_idx, , drop = FALSE] # next-month daily volatility

    # Hook for nonlinear shrinkage later
    # Long run target S
    S_target <- if (is.null(S_builder)) {
      safe_cor(Z_insample)
    } else {
      S_builder(Z_insample)
    }

    # Estimate DCC model on in-sample block, if it fails just retunr  null instead
    # of crashing everything
    dcc_fit <- tryCatch(
      estimate_dcc(Z_insample, S_target = S_target),
      error = function(e) NULL
    )

    if (is.null(dcc_fit)) {
      out[[i]] <- list(
        refit_date = blk$refit_date,
        forecast_dates = blk$forecast_dates,
        a = NA_real_,
        b = NA_real_,
        R_forecasts = NULL,
        H_forecasts = NULL,
        H_month = NULL
      )
      next
    }

    # forecast the full sequence of next-month daily conditional correlations
    dcc_fc <- forecast_dcc_correlations(
      Q_T = dcc_fit$Q_T,
      S = dcc_fit$S,
      a = dcc_fit$a,
      b = dcc_fit$b,
      h = length(blk$forecast_idx)
    )

    H_fc <- build_cov_from_sigma_and_R(
      SIGMA_future = SIGMA_future,
      R_list = dcc_fc,
      dates_future = blk$forecast_dates
    )

    # Main aggregated matrix used for precision matrix construction later
    H_month <- aggregate_monthly_cov(H_fc)

    # Stores everything from that monthly estimation/forecast cycle
    out[[i]] <- list(
      refit_month = blk$refit_month,
      refit_date = blk$refit_date,
      insample_idx = blk$insample_idx,
      forecast_idx = blk$forecast_idx,
      forecast_dates = blk$forecast_dates,
      a = dcc_fit$a,
      b = dcc_fit$b,
      S = dcc_fit$S,
      Q_T = dcc_fit$Q_T,
      R_forecasts = lapply(dcc_fc, `[[`, "R"),
      H_forecasts = H_fc,
      H_month = H_month,
      opt = dcc_fit$opt
    )
  }

  names(out) <- sapply(out, function(x) as.character(x$refit_date))
  out
}

# Running small DCC sample to check
nrow(Z_block)
nrow(SIGMA_block)
length(managed_portfolios$date)
class(Z_block)
class(SIGMA_block)

n_test <- min(
  1100,
  nrow(Z_block),
  nrow(SIGMA_block),
  length(managed_portfolios$date)
)

system.time({
  dcc_test <- run_dcc_monthly(
    Z_block = Z_block[1:n_test, , drop = FALSE],
    SIGMA_block = SIGMA_block[1:n_test, , drop = FALSE],
    dates_block = managed_portfolios$date[1:n_test],
    min_corr_window = 252,
    S_builder = build_nlshrink_target,
    trace = TRUE
  )
})

length(dcc_test)

summary_df <- data.frame(
  refit_date = as.Date(sapply(dcc_test, `[[`, "refit_date")),
  a = sapply(dcc_test, `[[`, "a"),
  b = sapply(dcc_test, `[[`, "b"),
  has_H_month = sapply(dcc_test, function(x) !is.null(x$H_month))
)

print(summary_df)
print(summary(summary_df$a))
print(summary(summary_df$b))
print(table(summary_df$has_H_month))

first_ok <- which(sapply(dcc_test, function(x) {
  !is.null(x$H_month) && all(is.finite(x$H_month))
}))[1]

res <- dcc_test[[first_ok]]

res$refit_date
res$a
res$b
dim(res$H_month)
range(eigen(res$H_month, symmetric = TRUE)$values)

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



## This is just for veryfying the performance of the benchmarks

# Calculating the Sharpe Ratios
sharpe_ratios_benchmarks <- tibble(
  Benchmark = c("EW Buy-and-Hold", "MVE"),
  SR = c(round(compute_SR(EW_returns_df$monthly_return),4)*sqrt(12), 
         round(compute_SR(MVE_returns_df$MVE_strategy_return),4)*sqrt(12))
)


# Creating data frame for letter regressions
benchmarks_returns <- EW_returns_df |> left_join(MVE_returns_df, by = "month") 
colnames(benchmarks_returns) <- c("month", "EW", "MVE")




# ====================================
# Generating figures 
# ====================================


# Cumulative wealth dataframe

cumulative_wealth_df <- benchmarks_returns |>
  mutate(
    EW = cumprod(1 + EW),
    MVE = cumprod(1 + MVE)
  )

axis <- par(lab = c(20, 8, 5))
plot(x = cumulative_wealth_df$month, y = cumulative_wealth_df$EW, xlab = "Date", 
     ylab = "Cumulative return", type = "l", col = "black", lwd = 1, lty = 1, ylim = c(0,60))
y_ticks <- pretty(cumulative_wealth_df$MVE)
abline(h = y_ticks, col = "grey85", lty = 1)
lines(x = cumulative_wealth_df$month, y = cumulative_wealth_df$EW, col = "black", lwd = 1)
lines(x = cumulative_wealth_df$month, y = cumulative_wealth_df$MVE , col = "blue", lty = 2)
legend("topleft", legend = c("BH", "MVE", "NET"), col = c("black", "blue", "red"), lty = c(1,2, 3), lwd = 2, bty = "n", cex = 0.8)


# Rolling Sharpe Ratio
rolling_SR_df <- data.frame(
  date = benchmarks_returns$month[12:dim(benchmarks_returns)[1]],
  SR_EW = NA_real_,
  SR_MVE = NA_real_
)

for(i in 12: dim(benchmarks_returns)[1]){
  rolling_SR_df$SR_EW[i-11] <- mean(benchmarks_returns$EW[(i-11):i])/sd(benchmarks_returns$EW[(i-11):i]) * sqrt(12)
  rolling_SR_df$SR_MVE[i-11] <- mean(benchmarks_returns$MVE[(i-11):i])/sd(benchmarks_returns$MVE[(i-11):i]) * sqrt(12)
}


plot(x = rolling_SR_df$date, y = rolling_SR_df$SR_EW, type = "l", xlab = "Date", ylab = "Rolling SR", lty = 1, col = "black")
y_ticks <- pretty( rolling_SR_df$SR_EW)
abline(h = y_ticks, col = "grey85", lty = 1)
lines(x = rolling_SR_df$date, y = rolling_SR_df$SR_EW, type = "l", col = "black", lty = 1)
lines(x = rolling_SR_df$date, y = rolling_SR_df$SR_MVE, col = "blue", lty = 2)
legend("topright", legend = c("BH", "MVE", "NET"), col = c("black", "blue", "red"), lty = c(1,2, 3), bty = "n", lwd = 2, cex = 0.8)

# Drowndown figure
drawdown_df <- cumulative_wealth_df |>
  arrange(month) |>
  mutate(
    EW_drawdown = EW / cummax(EW) - 1,
    MVE_drawdown = MVE / cummax(MVE) - 1
  )


plot(x = drawdown_df$month, y = drawdown_df$EW_drawdown, type = "l", ylim = c(-0.6,0), 
     xlab = "Date", ylab = "Drawdown")
y_ticks <- pretty(drawdown_df$EW_drawdown)
abline(h = y_ticks, col = "grey85", lty = 1)
lines(x = drawdown_df$month, y = drawdown_df$EW_drawdown, type = "l", col = "black", lty = 1)
lines(x = drawdown_df$month, y = drawdown_df$MVE_drawdown, type = "l", col = "blue", lty = 2)
legend("bottomright", legend = c("BH", "MVE", "NET"), col = c("black", "blue", "red"), lty = c(1,2, 3), bty = "n", lwd = 2, cex = 0.8)



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
perf_EW <- performance_summary(benchmarks_monthly$EW)
perf_MVE <- performance_summary(benchmarks_monthly$MVE)

#just to make a nice table with the different strategy
performance_table <- bind_rows(
  `EW Buy-and-Hold` = perf_EW,
  `MVE Vol-Managed` = perf_MVE,
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

sharpe_difference <- function(r1, r2, scale = 12) {
  compute_SR(r1, scale) - compute_SR(r2, scale)
}
sr_diff <- sharpe_difference(
  benchmarks_monthly$MVE,
  benchmarks_monthly$EW
)