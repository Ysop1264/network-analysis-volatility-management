#### Uncomment if packages not installed
# install.packages("tidyverse")
# install.packages("tidyfinance")
# install.packages("scales")
# install.packages("frenchdata")
# install.packages("dplyr")
# install.packages("moments")
# install.packages("sandwich")
# install.packages("rlang")
# install.packages("lmtest")
# install.packages("lubridate")
# install.packages("nlshrink")
# install.packages("rugarch")
# install.packages("xts")
# install.packages("zoo")
# install.packages("igraph")
# install.packages("PeerPerformance")
# install.packages("xdcclarge")
# install.packages("purrr")
# install.packages("tidyr")
# install.packages("zoo")
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
library(xdcclarge)
library(tidyr)
library(purrr)
library(zoo)

# DATA download and transformation
start_date <- as.Date("1971-01-01")
end_date   <- as.Date("2026-01-31")

# Code from tidy-finance website to download data
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
#to run for just the industry portfolios
# managed_portfolios <- factors_joined_excess |> select(-risk_free, -mkt_excess, -smb, -hml, -rmw, -cma, -mom)
# print(head(managed_portfolios))

#managed_portfolios <- managed_portfolios[,1:7]





#### Note that this script is copied from the main script and that it is adjusted for the univariate managed portfolio consideration
# Thus, the comments after the DCC-NL part might bu somewhat disconnected from the actual code, as I have deleted someparts and added 
# small changes for the univariate, the code tho produces the tables and a figure for the appendix section on univariate considertation



#' Computes the Sharpe Ratio of the returns
#'
#'
#' @param r Vector of returns
#' @return Sharpe Ratio
compute_SR <- function(r){
  SR = mean(r)/sd(r)
  return(SR)
}




# Fitting univariate GARCH(1,1)
# Function to fit GARCH(1,1) for one return series on the estimation window

#'
#' Function that estimates GARCH(1,1) for each series
#'
#' @param returns dataframe with dates as the first column, and managed portfolio returns as second
#' @param est_window estimation window length
#' @param distribution distribution assumption for GARCH(1,1) model
#' @param garch_order garch order
#' @param arma_order arma model for mean
#'
#' @return list with parameters, fit object, and estimated objects: mean, volatility, residuals, standardized residuals
#'
garch_estimate_single <- function(
    returns,
    est_window = 504,
    distribution = "std",
    garch_order = c(1,1),
    arma_order = c(0,0)
){
  
  # Making sure of the form of returns and dates
  dates = as.Date(returns[[1]])
  returns <- as.numeric(returns[[2]])
  
  # Scaling returns
  returns <- 100 * returns
  
  # Critical checks for estimating GARCH model
  if(length(returns) < est_window){
    stop("GARCH(1,1) error: estimation window is longer than whole series of returns!")
  }
  
  
  # Setting the specifications of the univariate GARCH(1,1) model
  spec_model <- ugarchspec(
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
  
  # Fitting the garch on the estimation window, keeping rest as out.sample
  fit <- ugarchfit(
    spec = spec_model,
    data = returns,
    out.sample = length(returns) - est_window,
    solver = "hybrid"
    
  )
  
  # Creating the output (saved parameters, fitted volatilities, residuals,...)
  parameters <- coef(fit)
  conditional_mean <- fitted(fit)
  conditional_volatility <- sigma(fit)
  residuals <- residuals(fit)
  residuals_st <- residuals(fit, standardize = TRUE)
  
  est_results <- data.frame(
    
    mu = c(conditional_mean),
    sigma = c(conditional_volatility),
    residuals = c(residuals),
    residuals_standard = c(residuals_st)
  )
  
  results_list <- list(
    parameters_model = parameters,
    fit_object = fit,
    est_result = est_results
  )
  
  return(results_list)
}

#This function takes the whole dataset and runs GARCH individually for each assets
#'
#'Estimates GARCH(1,1) model for each return series
#' @param returns_df dataframe with dates as the first column, and managed portfolio as all other columns
#' @param est_window estimation window length
#' @param distribution distribution assumption for GARCH(1,1) model
#'
#' @return list of univariate garch(1,1) results for each series
#'
estimate_all_univariate_garch <- function(returns_df,
                                          est_window = 504,
                                          distribution = "norm"
){
  results <- lapply(colnames(returns_df)[-1], function(col_name) {
    
    # Create individual series data frame
    individual_series_df <- data.frame(
      date = as.Date(returns_df[[1]]),
      returns = returns_df[[col_name]]
    )
    
    garch_estimate_single(individual_series_df, est_window = est_window, distribution = distribution)
  })
  
  names(results) <- colnames(returns_df)[-1]
  results
}


#'
#'
#'Function that creates the inputs for DCC estimation from univariate garch models
#'
#' @param returns_df dataframe with dates as the first column, and managed portfolio as all other columns
#' @param univariate_garch_models list with all results from univariate garch models
#'
#' @return list contianing the H and Z matrix
#'
creating_inputs_for_DCC <-  function(returns_df, univariate_garch_models){
  H <- do.call(cbind,
               lapply(univariate_garch_models, function(x) x$est_result$sigma^2))
  colnames(H) <- colnames(returns_df)[-1]
  
  residuals_std <- do.call(cbind,
                           lapply(univariate_garch_models, function(x) x$est_result$residuals_standard))
  
  residuals_raw <- do.call(cbind,
                           lapply(univariate_garch_models, function(x) x$est_result$residuals))
  colnames(residuals_std) <- colnames(returns_df)[-1]
  colnames(residuals_raw) <- colnames(returns_df)[-1]
  
  list(
    H_initial = H,
    residuals_raw = residuals_raw,
    residuals_std = residuals_std
  )
  
}

#'
#'Function that creates the target matrix with non-linear shrinkage method
#'This produces C_hat, the NL-shrunk unconditional correlation target
#'for the DCC dynamics (Eq. 3 in proposal)
#'
#'@param residuals_matrix matrix of standardized returns for DCC estimation
#'
#'@return C target matrix
create_target_matrix <- function(residuals_matrix){
  S <- nlshrink_cov(residuals_matrix)
  
  # if (!is.matrix(S) || nrow(S) != ncol(S)) {
  #   stop("nlshrink_cov did not return a square matrix")
  # }
  
  C <- cov2cor(S)
  C <- make_psd(C)
  diag(C) <- 1
  
  return(C)
}


# Function that estimates DCC parameters using the xdcclarge package
# Falls back to manual optimisation if the package fails
# @param dcc_inputs output from creating_inputs_for_DCC()
# @param target_C the NL-shrinkage target (used as fallback)
# @return list with alpha, beta, and the target matrix used
estimate_dcc_parameters <- function(dcc_inputs, target_C = NULL) {
  
  H_mat <- as.matrix(dcc_inputs$H_initial)
  Z_mat <- as.matrix(dcc_inputs$residuals_raw)
  
  # Standardized residuals for manual DCC fallback
  Z_std <- Z_mat / sqrt(H_mat)
  Z_std[!is.finite(Z_std)] <- 0
  
  # Try using xdcclarge with NL shrinkage method first
  dcc_fit <- tryCatch({
    cdcc_estimation(
      ini.para = c(0.05, 0.93),
      ht = H_mat,
      residuals = Z_mat,
      method = "NLS"
    )
  }, error = function(e) {
    message("xdcclarge NLS method failed: ", e$message)
    message("Falling back to manual DCC estimation with NL target...")
    NULL
  })
  
  if (!is.null(dcc_fit)) {
    print(class(dcc_fit))
    print(names(dcc_fit))
    str(dcc_fit, max.level = 2)
  }
  
  if (!is.null(dcc_fit)) {
    alpha_hat <- NULL
    beta_hat  <- NULL
    
    # Direct slots, if they exist
    if (!is.null(dcc_fit$a) && !is.null(dcc_fit$b) &&
        length(dcc_fit$a) == 1 && length(dcc_fit$b) == 1 &&
        is.finite(dcc_fit$a) && is.finite(dcc_fit$b)) {
      alpha_hat <- as.numeric(dcc_fit$a)
      beta_hat  <- as.numeric(dcc_fit$b)
    }
    
    # Top-level parameters, if it exists
    if ((is.null(alpha_hat) || is.null(beta_hat)) &&
        !is.null(dcc_fit$par) && length(dcc_fit$par) >= 2) {
      alpha_hat <- as.numeric(dcc_fit$par[1])
      beta_hat  <- as.numeric(dcc_fit$par[2])
    }    
    
    # result$par (list of parameters) , which is what your object actually uses
    if ((is.null(alpha_hat) || is.null(beta_hat)) &&
        !is.null(dcc_fit$result) &&
        !is.null(dcc_fit$result$par) &&
        length(dcc_fit$result$par) >= 2) {
      alpha_hat <- as.numeric(dcc_fit$result$par[1])
      beta_hat  <- as.numeric(dcc_fit$result$par[2])
    }    
    
    # parameter, just in case
    if ((is.null(alpha_hat) || is.null(beta_hat)) &&
        !is.null(dcc_fit$para) && length(dcc_fit$para) >= 2) {
      alpha_hat <- as.numeric(dcc_fit$para[1])
      beta_hat  <- as.numeric(dcc_fit$para[2])
    }
    
    if (is.null(alpha_hat) || is.null(beta_hat) ||
        length(alpha_hat) != 1 || length(beta_hat) != 1 ||
        !is.finite(alpha_hat) || !is.finite(beta_hat)) {
      message("xdcclarge returned object, but could not extract scalar alpha/beta.")
      message("Falling back to manual DCC estimation with NL target...")
      dcc_fit <- NULL
    }
    
    if (!is.null(dcc_fit)) {
      if (!is.null(target_C)) {
        target_used <- target_C
      } else {
        target_used <- create_target_matrix(Z_std)
      }
      
      return(list(
        alpha = as.numeric(alpha_hat),
        beta  = as.numeric(beta_hat),
        target_C = as.matrix(target_used),
        method = "xdcclarge_NLS"
      ))
    }
  }
  
  # If xdcclarge failed, fall back to manual optimisation with our own NL target
  # This is the same approach as the old code but simplified
  if (is.null(target_C)) {
    target_C <- create_target_matrix(Z_std)
  }
  
  # Negative log-likelihood for DCC (same as old code)
  dcc_negloglik <- function(par, Z, S) {
    a <- par[1]; b <- par[2]
    if (a < 0 || b < 0 || (a + b) >= 0.999) return(1e12)
    
    Tn <- nrow(Z)
    N <- ncol(Z)
    Q_prev <- S
    ll <- 0
    #updates the DCC correlation - driving matrix Q_t
    for (t in seq_len(Tn)) {
      zlag <- if (t == 1) rep(0, N) else Z[t-1, ]
      zlag[!is.finite(zlag)] <- 0
      
      Q_t <- (1 - a - b) * S + a * tcrossprod(zlag) + b * Q_prev
      Q_t <- (Q_t + t(Q_t)) / 2
      
      # Normalising Q to get R which is the conditional correlation matrix 
      qd <- sqrt(pmax(diag(Q_t), 1e-10))
      Dinv <- diag(1/qd, nrow = length(qd))
      R_t <- Dinv %*% Q_t %*% Dinv
      diag(R_t) <- 1
      
      # Ridge regularisation
      R_reg <- R_t + diag(1e-8, N) 
      
      # Cholesky decomposition for faster computation
      U <- tryCatch(chol(R_reg), error = function(e) NULL)
      if (is.null(U)) return(1e12)
      
      logdet <- 2 * sum(log(diag(U)))
      if (!is.finite(logdet)) return(1e12)
      
      zt <- Z[t, ]
      if (any(!is.finite(zt))) return(1e12)
      
      y <- tryCatch(backsolve(U, zt, transpose = TRUE), error = function(e) NULL)
      if (is.null(y) || any(!is.finite(y))) return(1e12)
      
      x <- tryCatch(backsolve(U, y), error = function(e) NULL)
      if (is.null(x) || any(!is.finite(x))) return(1e12)
      
      ll <- ll + logdet + sum(zt * x)
      Q_prev <- Q_t
    }
    
    0.5 * ll
  }
  
  opt <- optim(
    par = c(0.05, 0.93),
    fn = function(p) dcc_negloglik(p, Z_std, target_C),
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    upper = c(0.5, 0.999)
  )
  
  if (opt$convergence != 0) {
    stop("Manual DCC optimization failed to converge.")
  }
  
  list(
    alpha = opt$par[1],
    beta  = opt$par[2],
    target_C = target_C,
    method = "manual_fallback"
  )
}


#'
#' Function that estimates the path of the Q matrix in-sample period
#' Needed to recover the terminal Q_T for forecasting
#'
#' @param alpha alpha from DCC estimate
#' @param beta beta from DCC estimate
#' @param target_C target unconditional correlation matrix
#' @param residuals_dcc matrix with in-sample residuals
#' @param est_window length of estimation window
#'
#' @return last estimated Q
#'
dcc_path <- function(alpha, beta, target_C, residuals_dcc, est_window){
  alpha <- as.numeric(alpha)
  beta  <- as.numeric(beta)
  target_C <- as.matrix(target_C)
  residuals_dcc <- as.matrix(residuals_dcc)
  
  if (length(alpha) != 1 || !is.finite(alpha)) stop("alpha must be one finite scalar")
  if (length(beta)  != 1 || !is.finite(beta)) stop("beta must be one finite scalar")
  if (!is.matrix(target_C) || nrow(target_C) != ncol(target_C)) stop("target_C must be square")
  if (ncol(residuals_dcc) != nrow(target_C)) stop("residuals_dcc columns must match target_C dimension")
  
  Q_prev <- target_C
  
  if (est_window <= 1) return(Q_prev)
  
  for(i in 2:est_window){
    zlag <- as.numeric(residuals_dcc[i - 1, ])
    zlag[!is.finite(zlag)] <- 0
    
    zz <- tcrossprod(zlag)
    Q_prev <- (1 - alpha - beta) * target_C + alpha * zz + beta * Q_prev
    Q_prev <- (Q_prev + t(Q_prev)) / 2
  }
  
  return(Q_prev)
}


# Function that aggregates the daily covariance forecasts to monthly covariance forecast
# Implements Eq. 6 in the proposal
# @param covariance_matricies list of daily forecasted matricies in the month
dcc_aggregate_daily_to_monthly <- function(covariance_matricies){
  monthly_covariance_matrix <- as.matrix(Reduce('+', covariance_matricies))
}


# Function that forecasts the monthly covariance matrix from the DCC estimates
# For h=1: uses last in-sample residual
# For h>1: uses mean-reversion formula Q_h = (1-a-b)*C + (a+b)*Q_{h-1}
# Then reconstructs H = D * R * D and aggregates to monthly
# @param alpha estimated alpha for DCC
# @param beta estimated beta for DCC
# @param final_Q last estimated Q from the in-sample period
# @param target_C target unconditional correlation matrix from the in-sample
# @param residuals_dcc matrix with in-sample residuals
# @param no_days_month number of trading days in the considered month
# @param D_t matrix of forecasted volatilities from univariate GARCH
# @return estimated monthly covariance matrix
dcc_forecast_monthly <- function(alpha, beta, final_Q, target_C, residuals_dcc, no_days_month, D_t){
  
  output_Q <- list()
  output_R <- list()
  output_H <- list()
  
  for(i in 1:no_days_month){
    
    if(i == 1){
      output_Q[[i]] <- (1 - alpha - beta) * target_C +
        alpha * tcrossprod(as.numeric(residuals_dcc[nrow(residuals_dcc), ])) + beta * final_Q
      output_R[[i]] <- cov2cor(output_Q[[i]])
      output_H[[i]] <- D_t[[i]] %*% output_R[[i]] %*% D_t[[i]]
    } else {
      output_Q[[i]] <- (1 - alpha - beta) * target_C +
        (alpha + beta) * output_Q[[i-1]]
      output_R[[i]] <- cov2cor(output_Q[[i]])
      output_H[[i]] <- D_t[[i]] %*% output_R[[i]] %*% D_t[[i]]
    }
  }
  
  forecasted_monthly_cov_mat <- dcc_aggregate_daily_to_monthly(output_H)
  return(forecasted_monthly_cov_mat)
}


# Forecasts univariate garch volatilities using the already-fitted model object
# Avoids having to refit the GARCH again
# @param fit_object the ugarchfit object from garch_estimate_single()
# @param fc_window number of days ahead to forecast
# @returns data frame with mu and sigma forecasts
forecast_from_fitted_garch <- function(fit_object, fc_window = 21) {
  
  fc <- ugarchforecast(
    fitORspec = fit_object,
    n.ahead = fc_window,
    n.roll = 0
  )
  
  fc_results <- data.frame(
    mu = c(fitted(fc)),
    sigma = c(sigma(fc))
  )
  
  return(fc_results)
}

# Function that estimates the diagonal matrix of conditional volatilities based on univariate garch models
# Uses the already fitted GARCH objects instead of refitting
# @param univariate_garch_models list from estimate_all_univariate_garch()
# @param fc_window forecasting window (number of days in next month)
# @return list of D_t matricies for forecasted period length
estimate_D_t <- function(univariate_garch_models,
                         fc_window = 21
){
  garch_forecasts <- lapply(univariate_garch_models, function(x){
    forecast_from_fitted_garch(x$fit_object, fc_window = fc_window)
  })
  
  
  sigma_mat <- do.call(
    cbind,
    lapply(garch_forecasts, function(x) x$sigma)
  )
  
  colnames(sigma_mat) <- names(univariate_garch_models)
  
  D_list <- lapply(seq_len(nrow(sigma_mat)), function(t) {
    diag(sigma_mat[t, ])
  })
  
  return(list(
    sigma_mat = sigma_mat,
    D_list = D_list
  ))
  
}


# Forcing matrix to be positive semi-definite
# Done to allow inversion, valid likelihood and prevent DCC explosion
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

# Function to convert covariance matrix to correlation
cov_to_cor <- function(Sigma, eps = 1e-10) {
  d <- sqrt(pmax(diag(Sigma), eps))
  Dinv <- diag(1 / d, nrow = length(d))
  Dinv %*% Sigma %*% Dinv
}

# Adjacency matrix construction
# @param sigma_hat is the forecasted covariance matrix
# @return List with partial correlations and the adjacency matrix
adjacency_matrix <- function(sigma_hat, tau = 0.05, ridge = 1e-6,
                             use_correlation = TRUE) {
  sigma_hat <- as.matrix(sigma_hat)
  # Check to prevent singularity
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








# Monthly Index Helpers
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




# Two functions for creaing the volatility and volatility and centrality scalling (formulas in the appendix)
# It is the univariate testing case, univariate in the sense that we test how incorporating centrality for volatility timing for one factor/managed portfolio imporves (or not) the performance
compute_centrality_scale <- function(sigma_vec = NULL, ec_pos = NULL, ec_neg = NULL, lambda_pos = 0.1, lambda_neg = 0.1, eps = 1e-4){
  c = 1 + lambda_pos * ec_pos - lambda_neg * ec_neg
  return(1/(pmax(c, eps) * (sigma_vec^2)))
}

compute_var_scale <- function(sigma_vec = NULL){
  return(1/ (sigma_vec^2))
}


# Master function
# Loops over month-end refit dates and does the full monthly exercise:
# 1. Estimate univariate GARCH on in-sample
# 2. Build Z and H for DCC
# 3. Compute NL-shrinkage target
# 4. Estimate DCC (alpha, beta) via xdcclarge
# 5. Run DCC path to get Q_T
# 6. Forecast D_t for next month
# 7. Forecast DCC correlations -> monthly covariance
# 8. Build network: adjacency, centrality, spillovers, weights
# 9. Compute portfolio return for next month
run_dcc_network_monthly <- function(managed_portfolios,
                                    est_window = 504,
                                    min_corr_window = 252,
                                    rolling_window = NULL,
                                    distribution = "norm",
                                    should_reestimate = FALSE,
                                    reestimation_period = 12,
                                    tau = 0.05,
                                    lambda_pos = 0.1,
                                    lambda_neg = 0.1,
                                    eps = 1e-4,
                                    trace = TRUE) {
  
  dates <- as.Date(managed_portfolios$date)
  asset_names <- colnames(managed_portfolios)[-1]
  N <- length(asset_names)
  
  # Getting the monthly rebalancing schedule
  month_map <- make_month_index(
    dates = dates,
    min_obs = max(est_window, min_corr_window),
    rolling_window = rolling_window
  )
  
  if (length(month_map) == 0) {
    stop("No eligible monthly refit dates found after est_window.")
  }
  
  # Storage
  out <- vector("list", length(month_map))
  
  # Frozen parameters for re-estimation control
  frozen_alpha <- NULL
  frozen_beta <- NULL
  frozen_target_C <- NULL
  frozen_garch_models <- NULL
  months_since_estimation <- 0
  
  for (i in seq_along(month_map)) {
    
    blk <- month_map[[i]]
    fc_window <- length(blk$forecast_idx)
    
    if (trace) {
      message(
        "Month ", i, "/", length(month_map),
        " | Refit date: ", as.character(blk$refit_date),
        " | insample n = ", length(blk$insample_idx),
        " | forecast h = ", fc_window
      )
    }
    
    # Deciding whether to re-estimate or use frozen parameters
    need_estimation <- is.null(frozen_alpha) ||
      (should_reestimate && months_since_estimation >= reestimation_period)
    
    if (need_estimation) {
      
      if (trace) message("  >> Estimating GARCH + DCC parameters...")
      
      # In-sample data slice
      insample_df <- managed_portfolios[blk$insample_idx, ]
      
      # Step 1-2: Fit univariate GARCH(1,1) for all series
      garch_models <- tryCatch(
        estimate_all_univariate_garch(insample_df, est_window = est_window,
                                      distribution = distribution),
        error = function(e) {
          message("  GARCH fitting failed at ", as.character(blk$refit_date), ": ", e$message)
          NULL
        }
      )
      
      if (is.null(garch_models)) {
        out[[i]] <- list(refit_date = blk$refit_date, H_month = NULL)
        next
      }
      
      # Step 3: Build DCC inputs
      dcc_inputs <- creating_inputs_for_DCC(insample_df, garch_models)
      
      # Step 4: NL-shrinkage target
      target_C <- tryCatch(
        create_target_matrix(dcc_inputs$residuals_std),
        error = function(e) {
          message("  NL-shrinkage target failed at ", as.character(blk$refit_date), ": ", e$message)
          NULL
        }
      )
      
      print("TARGET C DIM:")
      print(dim(target_C))
      print("TARGET C SAMPLE:")
      print(target_C[1:3, 1:3])
      
      if (is.null(target_C)) {
        out[[i]] <- list(refit_date = blk$refit_date, H_month = NULL)
        next
      }
      
      # Step 5: Estimate DCC parameters via xdcclarge (with fallback)
      dcc_est <- tryCatch(
        estimate_dcc_parameters(dcc_inputs, target_C),
        error = function(e) {
          message("  DCC estimation failed at ", as.character(blk$refit_date), ": ", e$message)
          NULL
        }
      )
      
      print(dcc_est$alpha)
      print(dcc_est$beta)
      print(dcc_est$method)
      
      if (is.null(dcc_est)) {
        out[[i]] <- list(refit_date = blk$refit_date, H_month = NULL)
        next
      }
      
      # Freeze the parameters
      frozen_alpha <- dcc_est$alpha
      frozen_beta  <- dcc_est$beta
      frozen_target_C <- dcc_est$target_C
      frozen_garch_models <- garch_models
      months_since_estimation <- 0
      
      if (trace) {
        message(sprintf("  >> DCC params: alpha=%.6f, beta=%.6f (method=%s)",
                        frozen_alpha, frozen_beta, dcc_est$method))
      }
      
    } else {
      
      # Use the last jointly estimated GARCH + DCC objects until next scheduled re-estimation
      garch_models <- frozen_garch_models
      
      if (is.null(garch_models)) {
        message("  No frozen GARCH models available.")
        out[[i]] <- list(refit_date = blk$refit_date, H_month = NULL)
        next
      }
      
      # Keep using the in-sample slice only for bookkeeping if needed
      insample_df <- managed_portfolios[blk$insample_idx, ]
      
      # Rebuild DCC inputs from the frozen GARCH results
      dcc_inputs <- creating_inputs_for_DCC(insample_df, garch_models)
      
      months_since_estimation <- months_since_estimation + 1
    }
    
    # Step 6: Run DCC path to recover terminal Q_T
    Z_insample <- dcc_inputs$residuals_std
    
    Q_T <- dcc_path(
      alpha = frozen_alpha,
      beta = frozen_beta,
      target_C = frozen_target_C,
      residuals_dcc = Z_insample,
      est_window = nrow(Z_insample)
    )
    
    # Step 7: Forecast D_t for the next month using fitted GARCH objects
    D_t_result <- tryCatch(
      estimate_D_t(garch_models, fc_window = fc_window),
      error = function(e) {
        message("  D_t forecast failed at ", as.character(blk$refit_date), ": ", e$message)
        NULL
      }
    )
    
    if (is.null(D_t_result)) {
      out[[i]] <- list(refit_date = blk$refit_date, H_month = NULL)
      next
    }
    print(dim(Q_T))
    print(dim(frozen_target_C))
    print(length(asset_names))
    # Step 8: Forecast monthly covariance matrix
    H_month <- tryCatch(
      dcc_forecast_monthly(
        alpha = frozen_alpha,
        beta = frozen_beta,
        final_Q = Q_T,
        target_C = frozen_target_C,
        residuals_dcc = Z_insample,
        no_days_month = fc_window,
        D_t = D_t_result$D_list
      ),
      error = function(e) {
        message("  DCC forecast failed at ", as.character(blk$refit_date), ": ", e$message)
        NULL
      }
    )
    
    if (is.null(H_month)) {
      out[[i]] <- list(refit_date = blk$refit_date, H_month = NULL)
      next
    }
    
    colnames(H_month) <- asset_names
    rownames(H_month) <- asset_names
    
    # Step 9: Network construction
    sigma_vec <- sqrt(pmax(diag(H_month), 0))
    print(summary(sigma_vec))
    names(sigma_vec) <- asset_names
    
    adj <- adjacency_matrix(
      sigma_hat = H_month,
      tau = tau,
      ridge = 1e-6,
      use_correlation = TRUE
    )
    
    cent <- network_centrality(
      adjacency_pos = adj$adjacency_pos,
      adjacency_neg = adj$adjacency_neg
    )
    
 
    
    # Step 10: Compute realized return for next month
    next_month_returns <- managed_portfolios[blk$forecast_idx, -1, drop = FALSE]
    realized_monthly <- apply(next_month_returns, 2, function(z) prod(1 + z, na.rm = TRUE) - 1)
    
    
    # the univariate strategies 
    var_scale <- compute_var_scale(sigma_vec = sigma_vec)
    centrality_scale <-compute_centrality_scale(sigma_vec = sigma_vec, ec_pos = cent$ec_pos, ec_neg = cent$ec_neg, lambda_pos = lambda_pos, lambda_neg = lambda_neg)

    univ_var_return <- var_scale * realized_monthly
    univ_centrality <- centrality_scale * realized_monthly
    

    
    
    # Store everything
    out[[i]] <- list(
      refit_month = blk$refit_month,
      refit_date = blk$refit_date,
      insample_idx = blk$insample_idx,
      forecast_idx = blk$forecast_idx,
      forecast_dates = blk$forecast_dates,
      alpha = frozen_alpha,
      beta = frozen_beta,
      H_month = H_month,
      sigma_vec = setNames(as.numeric(sigma_vec), asset_names),
      partial_corr = `dimnames<-`(adj$partial_corr, list(asset_names, asset_names)),
      adjacency_pos = `dimnames<-`(adj$adjacency_pos, list(asset_names, asset_names)),
      adjacency_neg = `dimnames<-`(adj$adjacency_neg, list(asset_names, asset_names)),
      ec_pos = setNames(as.numeric(cent$ec_pos), asset_names),
      ec_neg = setNames(as.numeric(cent$ec_neg), asset_names),
      var_scale = setNames(as.numeric(var_scale), asset_names),
      centrality_scale = setNames(as.numeric(centrality_scale), asset_names),
      univ_var_return = setNames(as.numeric(univ_var_return), asset_names),
      univ_centrality = setNames(as.numeric(univ_centrality), asset_names),
      realized_returns = setNames(as.numeric(realized_monthly), asset_names)
    )
  }
  
  
  
  
  
  
  names(out) <- sapply(out, function(x) as.character(x$refit_date))
  
 
  # Keep only successful monthly outputs (should be all the months)
  ok <- sapply(out, function(z) !is.null(z$univ_var_return))
  
  
  # Double check and inform the terminal
  if (!any(ok)) stop("No valid monthly outputs produced.")
  
  
  dates_out <- as.Date(names(out)[ok])
  
  # Build factor-by-factor matrices
  sigma_forecast_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$sigma_vec)))
  colnames(sigma_forecast_mat) <- asset_names
  rownames(sigma_forecast_mat) <- names(out)[ok]
  
  ec_pos_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$ec_pos)))
  colnames(ec_pos_mat) <- asset_names
  rownames(ec_pos_mat) <- names(out)[ok]
  
  ec_neg_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$ec_neg)))
  colnames(ec_neg_mat) <- asset_names
  rownames(ec_neg_mat) <- names(out)[ok]
  
  var_scale_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$var_scale)))
  colnames(var_scale_mat) <- asset_names
  rownames(var_scale_mat) <- names(out)[ok]
  
  centrality_scale_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$centrality_scale)))
  colnames(centrality_scale_mat) <- asset_names
  rownames(centrality_scale_mat) <- names(out)[ok]
  
  realized_returns_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$realized_returns)))
  colnames(realized_returns_mat) <- asset_names
  rownames(realized_returns_mat) <- names(out)[ok]
  
  univ_var_return_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$univ_var_return)))
  colnames(univ_var_return_mat) <- asset_names
  rownames(univ_var_return_mat) <- names(out)[ok]
  
  univ_centrality_return_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$univ_centrality)))
  colnames(univ_centrality_return_mat) <- asset_names
  rownames(univ_centrality_return_mat) <- names(out)[ok]
  
  # Building the dataframes instead, for the later analusis (table and figures)
  sigma_forecast_df <- data.frame(
    date = dates_out,
    sigma_forecast_mat,
    
    
    # Using the names of the individual managed portfolios instead
    check.names = FALSE
  )
  
  ec_pos_df <- data.frame(
    date = dates_out,
    ec_pos_mat,
    check.names = FALSE
  )
  
  ec_neg_df <- data.frame(
    date = dates_out,
    ec_neg_mat,
    check.names = FALSE
  )
  
  
  var_scale_df <- data.frame(
    date = dates_out,
    var_scale_mat,
    check.names = FALSE
  )
  
  centrality_scale_df <- data.frame(
    date = dates_out,
    centrality_scale_mat,
    check.names = FALSE
  )
  
  realized_returns_df <- data.frame(
    date = dates_out,
    realized_returns_mat,
    check.names = FALSE
  )
  
  
  
  
  univ_var_return_df <- data.frame(
    date = dates_out,
    univ_var_return_mat,
    check.names = FALSE
  )
  
  univ_centrality_return_df <- data.frame(
    date = dates_out,
    univ_centrality_return_mat,
    check.names = FALSE
  )
  
  summary_df <- data.frame(
    period = seq_along(dates_out),
    date = dates_out,
    stringsAsFactors = FALSE
  )
  
  return(list(
    summary = summary_df,
    sigma_forecast = sigma_forecast_mat,
    ec_pos = ec_pos_mat,
    ec_neg = ec_neg_mat,
    var_scale = var_scale_mat,
    centrality_scale = centrality_scale_mat,
    realized_returns = realized_returns_mat,
    univariate_var_timed = univ_var_return_mat,
    univariate_var_centrality_timed = univ_centrality_return_mat,
    sigma_forecast_df = sigma_forecast_df,
    ec_pos_df = ec_pos_df,
    ec_neg_df = ec_neg_df,
    var_scale_df = var_scale_df,
    centrality_scale_df = centrality_scale_df,
    realized_returns_df = realized_returns_df,
    univ_var_return_df = univ_var_return_df,
    univ_centrality_return_df = univ_centrality_return_df,
    per_period = out
  ))
}

  





time_full_pipeline <- system.time({
  net_results_full <- run_dcc_network_monthly(
    managed_portfolios = managed_portfolios,
    est_window = 504,
    min_corr_window = 252,
    rolling_window = 504,
    distribution = "std",
    should_reestimate = TRUE,
    reestimation_period = 12, #24 # Rerunning for annual
    tau = sqrt(0.025),
    lambda_pos = 0.1,
    lambda_neg = 0.1,
    eps = 1e-4,
    trace = TRUE
  )
})

print(time_full_pipeline)










# Results before the scalling 
bh_mat_raw   <- as.matrix(net_results_full$realized_returns)
var_mat_raw  <- as.matrix(net_results_full$univariate_var_timed)
cent_mat_raw <- as.matrix(net_results_full$univariate_var_centrality_timed)

asset_names <- colnames(bh_mat_raw)

# Use the end date of the forecast month, not the refit month for the visualization and table
dates_monthly <- as.Date(sapply(net_results_full$per_period, function(z) {
  max(as.Date(z$forecast_dates))
}))



#Normalize strategy to match buy-and-hold volatility (the managed portfolio specifc c as in the appendix)
normalize_to_target_vol <- function(timed_mat, target_mat) {
  
  
  timed_mat  <- as.matrix(timed_mat)
  target_mat <- as.matrix(target_mat)
  
  c_vec <- sapply(seq_len(ncol(timed_mat)), function(j) {
    
    # added the na.rm but it should not contain any missing values (one we have non missing data and two there are already additional checks perfomred earlier)
    target_sd <- sd(target_mat[, j], na.rm = TRUE)
    timed_sd  <- sd(timed_mat[, j], na.rm = TRUE)
    
    # The managed portfolio specific c_i
    target_sd / timed_sd
  })
  
  scaled_mat <- sweep(timed_mat, 2, c_vec, `*`)
  colnames(scaled_mat) <- colnames(timed_mat)
  
  list(
    scaled = scaled_mat,
    c_vec = c_vec
  )
}


# Scaling the strategies
var_norm  <- normalize_to_target_vol(var_mat_raw,  bh_mat_raw)
cent_norm <- normalize_to_target_vol(cent_mat_raw, bh_mat_raw)


# Getting them out 
bh_mat   <- bh_mat_raw
var_mat  <- var_norm$scaled
cent_mat <- cent_norm$scaled



# Redoing the dataframe
bh_df <- data.frame(date = dates_monthly, bh_mat, check.names = FALSE)
var_df <- data.frame(date = dates_monthly, var_mat, check.names = FALSE)
cent_df <- data.frame(date = dates_monthly, cent_mat, check.names = FALSE)



# Annualized and percentages statistics for the table 
ann_mean <- function(x) {
  mean(x, na.rm = TRUE) * 12 * 100
}
ann_vol  <- function(x) {
  sd(x, na.rm = TRUE) * sqrt(12) * 100
}

ann_sr   <- function(x) {
  s <- sd(x, na.rm = TRUE)
  mean(x, na.rm = TRUE) / s * sqrt(12)
}






alpha_stats <- function(y, x) {
  df <- data.frame(y = y, x = x)
  df <- df[complete.cases(df), ]
  
  fit <- lm(y ~ x, data = df)
  
  
  nw  <- coeftest(fit, vcov. = NeweyWest(fit, prewhite = FALSE, adjust = TRUE))
  
  tibble(
    alpha = unname(coef(fit)[1]) * 12 * 100,
    t_alpha = unname(nw[1, 3]),
    beta = unname(coef(fit)[2])
  )
}




# Computing the statistics for the managed portfolios (if you guys have idea we can add something else to compare but I think this should be decent)
factor_results <- map_dfr(asset_names, function(a) {
  bh_x   <- bh_mat[, a]
  var_x  <- var_mat[, a]
  cent_x <- cent_mat[, a]
  
  stats_var  <- alpha_stats(var_x, bh_x)
  stats_cent <- alpha_stats(cent_x, bh_x)
  
  alpha_var_val     <- stats_var$alpha[[1]]
  t_alpha_var_val   <- stats_var$t_alpha[[1]]
  
  alpha_cent_val    <- stats_cent$alpha[[1]]
  t_alpha_cent_val  <- stats_cent$t_alpha[[1]]
  
  tibble(
    asset = a,
    
    mean_var = ann_mean(var_x),
    vol_var  = ann_vol(var_x),
    sr_var   = ann_sr(var_x),
    alpha_var = alpha_var_val,
    t_alpha_var = t_alpha_var_val,
    
    mean_cent = ann_mean(cent_x),
    vol_cent  = ann_vol(cent_x),
    sr_cent   = ann_sr(cent_x),
    alpha_cent = alpha_cent_val,
    t_alpha_cent = t_alpha_cent_val,
    
    delta_mean  = ann_mean(cent_x) - ann_mean(var_x),
    delta_sr    = ann_sr(cent_x) - ann_sr(var_x),
    delta_alpha = alpha_cent_val - alpha_var_val
  )
})


# Function for each row of table 4 (similar to what we do for the newest table 1 of dataset summary - basically cross sectional stats. over the assets)
cs_summary_row <- function(metric_name, x) {
  tibble(
    Metric = metric_name,
    Mean   = mean(x, na.rm = TRUE),
    SD     = sd(x, na.rm = TRUE),
    P25    = quantile(x, 0.25, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    P75    = quantile(x, 0.75, na.rm = TRUE)
  )
}

# Just additional function to round the results, bc I forgot it earlier 
round_table <- function(df, digits = 2) {
  df %>%
    mutate(across(-Metric, ~ round(.x, digits)))
}


panel_A <- bind_rows(
  cs_summary_row("Annualised mean return", factor_results$mean_cent),
  cs_summary_row("Annualised volatility",  factor_results$vol_cent),
  cs_summary_row("Annualised Sharpe ratio", factor_results$sr_cent),
  cs_summary_row("Annualised alpha", factor_results$alpha_cent)
) %>%
  round_table(2)



panel_B <- bind_rows(
  cs_summary_row("Annualised mean return", factor_results$mean_var),
  cs_summary_row("Annualised volatility",  factor_results$vol_var),
  cs_summary_row("Annualised Sharpe ratio", factor_results$sr_var),
  cs_summary_row("Annualised alpha", factor_results$alpha_var)
) %>%
  round_table(2)



# Table 2 for top gainers and bottom gainers
top5 <- factor_results |>
  arrange(desc(delta_alpha))|>
  transmute(
    managed_portfolio = asset,
    alpha_vol = round(alpha_var, 2),
    alpha_centrality = round(alpha_cent, 2),
    delta_alpha = round(delta_alpha, 2),
    SR_vol = round(sr_var, 2),
    SR_centrality = round(sr_cent, 2),
    delta_SR = round(delta_sr, 2)
  ) |> slice(1:5)

bottom5 <- factor_results |>
  arrange(delta_alpha) |>
  transmute(
    managed_portfolio = asset,
    alpha_vol = round(alpha_var, 2),
    alpha_centrality = round(alpha_cent, 2),
    delta_alpha = round(delta_alpha, 2),
    SR_vol = round(sr_var, 2),
    SR_centrality = round(sr_cent, 2),
    delta_SR = round(delta_sr, 2)
  ) |> slice(1:5)






# Making the share of portfolios figure


# theme of ggplot
theme_report <- function() {
  theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
}



# Computingt he rolling Sharpe Ratio
rolling_sr <- function(x, window = 12) {
  zoo::rollapplyr(
    x,
    width = window,
    FUN = function(z) {
      s <- sd(z, na.rm = TRUE)
      mean(z, na.rm = TRUE) / s * sqrt(12)
    },
    fill = NA
  )
}




window_months <- 12



roll_sr_var <- sapply(seq_len(ncol(var_mat)), function(j) {
  rolling_sr(var_mat[, j], window = window_months)
})
colnames(roll_sr_var) <- asset_names

roll_sr_cent <- sapply(seq_len(ncol(cent_mat)), function(j) {
  rolling_sr(cent_mat[, j], window = window_months)
})
colnames(roll_sr_cent) <- asset_names

delta_roll_sr <- roll_sr_cent - roll_sr_var


# rolling stats for the fig
rolling_sr_summary <- tibble(
  date = dates_monthly,
  p25 = apply(delta_roll_sr, 1, quantile, probs = 0.25, na.rm = TRUE),
  median = apply(delta_roll_sr, 1, median, na.rm = TRUE),
  p75 = apply(delta_roll_sr, 1, quantile, probs = 0.75, na.rm = TRUE),
  mean = rowMeans(delta_roll_sr, na.rm = TRUE),
  
  # This is for the graph below
  share_positive = rowMeans(delta_roll_sr > 0, na.rm = TRUE)
)



# Actual plotting

# Just a note for whoever is running: I supressed the warnings bc we have 11 NA values - as I need the window length 
fig_share_positive <- ggplot(rolling_sr_summary, aes(x = date, y = share_positive)) +
  geom_line(linewidth = 0.7, colour = "black", na.rm = TRUE) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.6) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "",
    subtitle = "",
    x = NULL,
    y = "Share with positive gain"
  ) +
  theme_report()

print(fig_share_positive)

ggsave(
  filename = "fig_share_positive_wide.png",
  plot = fig_share_positive,
  width = 14,
  height = 5,
  dpi = 300
)

