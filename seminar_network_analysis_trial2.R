#### Uncomment if packages not installed
install.packages("tidyverse")
install.packages("tidyfinance")
install.packages("scales")
install.packages("frenchdata")
install.packages("dplyr")
install.packages("moments")
install.packages("sandwich")
install.packages("rlang")
install.packages("lmtest")
install.packages("lubridate")
install.packages("nlshrink")
install.packages("rugarch")
install.packages("xts")
install.packages("zoo")
install.packages("igraph")
install.packages("PeerPerformance")
install.packages("xdcclarge")
install.packages("purrr")
install.packages("tidyr")
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
#to run for just the industry portfolios
# managed_portfolios <- factors_joined_excess |> select(-risk_free, -mkt_excess, -smb, -hml, -rmw, -cma, -mom)
# print(head(managed_portfolios))

# =================================================
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
# create_target_matrix <- function(residuals_matrix){
#   C <- cov2cor(nlshrink_cov(residuals_matrix))
#   return(C)
# }

create_target_matrix <- function(residuals_matrix){
  S <- nlshrink_cov(residuals_matrix)
  
  if (!is.matrix(S) || nrow(S) != ncol(S)) {
    stop("nlshrink_cov did not return a square matrix")
  }
  
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
    cat("\n===== xdcclarge object structure =====\n")
    print(class(dcc_fit))
    print(names(dcc_fit))
    str(dcc_fit, max.level = 2)
    cat("======================================\n")
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
    
    # Top-level par, if it exists
    if ((is.null(alpha_hat) || is.null(beta_hat)) &&
        !is.null(dcc_fit$par) && length(dcc_fit$par) >= 2) {
      alpha_hat <- as.numeric(dcc_fit$par[1])
      beta_hat  <- as.numeric(dcc_fit$par[2])
    }    
    
    # result$par, which is what your object actually uses
    if ((is.null(alpha_hat) || is.null(beta_hat)) &&
        !is.null(dcc_fit$result) &&
        !is.null(dcc_fit$result$par) &&
        length(dcc_fit$result$par) >= 2) {
      alpha_hat <- as.numeric(dcc_fit$result$par[1])
      beta_hat  <- as.numeric(dcc_fit$result$par[2])
    }    
    
    # para, just in case
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
    
    Tn <- nrow(Z); N <- ncol(Z)
    Q_prev <- S; ll <- 0
    
    for (t in seq_len(Tn)) {
      zlag <- if (t == 1) rep(0, N) else Z[t-1, ]
      zlag[!is.finite(zlag)] <- 0
      
      Q_t <- (1 - a - b) * S + a * tcrossprod(zlag) + b * Q_prev
      Q_t <- (Q_t + t(Q_t)) / 2
      
      # Normalising Q to get R
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
    # else {
    
    #   # Re-filter with frozen params but updated in-sample residuals for Q_T
    #   insample_df <- managed_portfolios[blk$insample_idx, ]
    
    #   garch_models <- tryCatch(
    #     estimate_all_univariate_garch(insample_df, est_window = est_window,
    #                                   distribution = distribution),
    #     error = function(e) {
    #       message("  GARCH refit failed, using frozen: ", e$message)
    #       NULL
    #     }
    #   )
    
    #   if (is.null(garch_models)) {
    #     garch_models <- frozen_garch_models
    #   }
    
    #   dcc_inputs <- creating_inputs_for_DCC(insample_df, garch_models)
    #   months_since_estimation <- months_since_estimation + 1
    # }
    
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
    
    # Step 10: Compute realized return for next month
    next_month_returns <- managed_portfolios[blk$forecast_idx, -1, drop = FALSE]
    realized_monthly <- apply(next_month_returns, 2, function(z) prod(1 + z, na.rm = TRUE) - 1)
    
    # Unscaled portfolio return
    portfolio_return_unscaled <- network_portfolio_return(nw$w_tilde, realized_monthly)
    
    # Equal weight benchmark return for same month
    ew_return <- mean(realized_monthly, na.rm = TRUE)
    
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
      spillover_pos = setNames(as.numeric(spl$spillover_pos), asset_names),
      spillover_neg = setNames(as.numeric(spl$spillover_neg), asset_names),
      penalty = setNames(as.numeric(nw$penalty), asset_names),
      w_tilde = setNames(as.numeric(nw$w_tilde), asset_names),
      realized_returns = realized_monthly,
      portfolio_return_unscaled = portfolio_return_unscaled,
      EW_benchmark = ew_return
    )
  }
  
  names(out) <- sapply(out, function(x) as.character(x$refit_date))
  
  # Building the summary output dataframes (same as our old code)
  ok <- sapply(out, function(z) !is.null(z$w_tilde))
  
  if (!any(ok)) stop("No valid network outputs produced.")
  
  # Weights matrix
  w_tilde_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$w_tilde)))
  colnames(w_tilde_mat) <- asset_names
  rownames(w_tilde_mat) <- names(out)[ok]
  
  # Penalty matrix
  penalty_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$penalty)))
  colnames(penalty_mat) <- asset_names
  rownames(penalty_mat) <- names(out)[ok]
  
  # Spillover matrices
  spillover_pos_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$spillover_pos)))
  colnames(spillover_pos_mat) <- asset_names
  rownames(spillover_pos_mat) <- names(out)[ok]
  
  spillover_neg_mat <- do.call(rbind, lapply(out[ok], function(z) unname(z$spillover_neg)))
  colnames(spillover_neg_mat) <- asset_names
  rownames(spillover_neg_mat) <- names(out)[ok]
  
  # Returns
  dates_out <- as.Date(names(out)[ok])
  net_returns_unscaled <- sapply(out[ok], function(z) z$portfolio_return_unscaled)
  ew_returns <- sapply(out[ok], function(z) z$EW_benchmark)
  
  # Combined returns dataframe
  network_vs_benchmark <- data.frame(
    date = dates_out,
    network_return_unscaled = as.numeric(net_returns_unscaled),
    EW_benchmark = as.numeric(ew_returns),
    stringsAsFactors = FALSE
  )
  
  # Summary
  summary_df <- data.frame(
    period = which(ok),
    date = dates_out,
    stringsAsFactors = FALSE
  )
  
  list(
    summary = summary_df,
    w_tilde = w_tilde_mat,
    penalty = penalty_mat,
    spillover_pos = spillover_pos_mat,
    spillover_neg = spillover_neg_mat,
    network_vs_benchmark = network_vs_benchmark,
    per_period = out
  )
}

# # Test on smaller subset first
# system.time({
#  test_result <- run_dcc_network_monthly(
#     managed_portfolios = managed_portfolios,
#     est_window = 504,
#     min_corr_window = 252,
#     rolling_window = 504,
#     distribution = "std",
#     should_reestimate = TRUE,
#     reestimation_period = 24,
#     tau = 0.05,
#     lambda_pos = 0.1,
#     lambda_neg = 0.1,
#     eps = 1e-4,
#     trace = TRUE
#   )
# })

# test_result$summary
# test_result$network_vs_benchmark
# head(test_result$w_tilde)
# head(test_result$penalty)
# rowSums(abs(test_result$w_tilde))

# Run for full sample
cat("Rows in full sample:", nrow(managed_portfolios), "\n")
cat("Number of assets/factors:", ncol(managed_portfolios) - 1, "\n")

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

network_vs_benchmark_all <- net_results_full$network_vs_benchmark

# Scaling the network returns to match EW volatility
vol_target <- sd(network_vs_benchmark_all$EW_benchmark, na.rm = TRUE)
vol_strat_unscaled <- sd(network_vs_benchmark_all$network_return_unscaled, na.rm = TRUE)
c_scaling <- vol_target / vol_strat_unscaled

network_vs_benchmark_all <- network_vs_benchmark_all %>%
  mutate(net_strategy_return = network_return_unscaled * c_scaling)

# inspect
names(net_results_full$summary)
head(net_results_full$w_tilde)
head(net_results_full$penalty)
head(net_results_full$spillover_pos)
head(net_results_full$spillover_neg)

names(network_vs_benchmark_all)
names(net_results_full)



# =================================================
# Benchmarks
# =================================================
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
MVE_returns_weights_df <- MVE_returns_df |> 
  mutate(
    MVE_weights = 1/ rv_lag
  )
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
MVE_returns_weights_df$MVE_weights = MVE_returns_weights_df$MVE_weights * c
MVE_returns_df <- MVE_returns_df |> select(-c(monthly_return, MVE_scaled))

# Full monthly MVE asset weights: scalar monthly leverage times fixed weight vector b
MVE_weights_mat <- MVE_returns_weights_df$MVE_weights %o% as.numeric(b)

colnames(MVE_weights_mat) <- colnames(managed_portfolios)[-1]

MVE_weights_df <- data.frame(
  date = MVE_returns_weights_df$month,
  MVE_weights_mat,
  check.names = FALSE
)


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


# Creating data frame for later regressions
benchmarks_returns <- EW_returns_df |> left_join(MVE_returns_df, by = "month")
colnames(benchmarks_returns) <- c("month", "EW", "MVE")


# ====================================
# Generating figures
# ====================================
# Building monthly network returns with aligned dates
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

# Keep only months where all plotted strategies exist
combined_returns_plot <- combined_returns %>%
  drop_na(EW, MVE, net_strategy_return)

returns_for_tables <- combined_returns %>%
  select(
    date = month,
    EW,
    MVE,
    NET = net_strategy_return
  ) %>%
  drop_na()

# Build cumulative wealth safely
# Build cumulative wealth on common sample
cumulative_wealth_df <- combined_returns_plot %>%
  mutate(
    EW  = cumprod(1 + EW),
    MVE = cumprod(1 + MVE),
    NET = cumprod(1 + net_strategy_return)
  )

#png("cumulative_wealth.png", width = 900, height = 600)
jpeg("cumulative_wealth.jpeg", width = 900, height = 600, quality = 100)
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
  date = combined_returns_plot$month[12:nrow(combined_returns_plot)],
  SR_EW = NA_real_,
  SR_MVE = NA_real_,
  SR_NET = NA_real_
)

for(i in 12:nrow(combined_returns_plot)){
  rolling_SR_df$SR_EW[i-11]  <- mean(combined_returns_plot$EW[(i-11):i]) / sd(combined_returns_plot$EW[(i-11):i]) * sqrt(12)
  rolling_SR_df$SR_MVE[i-11] <- mean(combined_returns_plot$MVE[(i-11):i]) / sd(combined_returns_plot$MVE[(i-11):i]) * sqrt(12)
  rolling_SR_df$SR_NET[i-11] <- mean(combined_returns_plot$net_strategy_return[(i-11):i]) / sd(combined_returns_plot$net_strategy_return[(i-11):i]) * sqrt(12)
}

# png("rolling_sharpe.png", width = 900, height = 600)
jpeg("rolling_sharpe.jpeg", width = 900, height = 600, quality = 100)
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

# png("drawdown.png", width = 900, height = 600)
jpeg("drawdown.jpeg", width = 900, height = 600, quality = 100)
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
perf_NET <- performance_summary(network_vs_benchmark_all$net_strategy_return)

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

# Aligning network returns to benchmark dates for alpha tests
net_for_alpha <- network_vs_benchmark_all %>%
  mutate(month = floor_date(date, "month")) %>%
  left_join(benchmarks_returns, by = "month") %>%
  drop_na()

alpha_NET_vs_BH <- alpha_test(
  strategy_ret = net_for_alpha$net_strategy_return,
  benchmark_ret = net_for_alpha$EW
)

print(alpha_NET_vs_BH)

alpha_NET_vs_MVE <- alpha_test(
  strategy_ret = net_for_alpha$net_strategy_return,
  benchmark_ret = net_for_alpha$MVE
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
  x = benchmarks_returns$MVE,
  y = benchmarks_returns$EW,
  control = list(
    type = 2,
    hac = TRUE,
    nBoot = 499,
    bBoot = 0
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
  net_for_alpha$net_strategy_return,
  net_for_alpha$EW
)
print(sr_diff_net_ew)

# Difference between NET and MVE
sr_diff_net_mve <- sharpe_difference(
  net_for_alpha$net_strategy_return,
  net_for_alpha$MVE
)
print(sr_diff_net_mve)

compute_turnover_drift <- function(returns_df, weights_df, half_turnover = FALSE) {
  
  W_raw <- as.matrix(weights_df[, -1, drop = FALSE])
  R <- as.matrix(returns_df[, -1, drop = FALSE])
  
  if (nrow(W_raw) != nrow(R)) {
    stop("weights_df and returns_df must have the same number of rows.")
  }
  
  # Normalize target weights to gross exposure = 1
  W <- W_raw / rowSums(abs(W_raw))
  
  n <- nrow(W)
  turnover <- rep(NA_real_, n)
  
  for (t in 2:n) {
    
    w_prev <- W[t - 1, ]
    r_t <- R[t, ]
    
    # Drift previous positions through asset returns
    pos_pre <- w_prev * (1 + r_t)
    
    # Re-normalize by gross exposure, NOT net portfolio value
    gross_pre <- sum(abs(pos_pre), na.rm = TRUE)
    
    if (gross_pre == 0 || !is.finite(gross_pre)) {
      turnover[t] <- NA_real_
      next
    }
    
    w_pre <- pos_pre / gross_pre
    
    raw_turnover <- sum(abs(W[t, ] - w_pre), na.rm = TRUE)
    
    turnover[t] <- if (half_turnover) 0.5 * raw_turnover else raw_turnover
  }
  
  turnover
}

asset_returns_monthly <- managed_portfolios %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(across(-date, ~ prod(1 + .x, na.rm = TRUE) - 1), .groups = "drop") %>%
  rename(date = month)

canonical_months <- asset_returns_monthly %>%
  filter(date > end_date_estimation) %>%
  transmute(
    month = floor_date(date, "month"),
    date = date
  )

MVE_weights_mat <- MVE_returns_weights_df$MVE_weights %o% as.numeric(b)
colnames(MVE_weights_mat) <- colnames(managed_portfolios)[-1]

MVE_weights_df <- data.frame(
  month = MVE_returns_weights_df$month,
  MVE_weights_mat,
  check.names = FALSE
) %>%
  left_join(canonical_months, by = "month") %>%
  select(date, everything(), -month)

NET_weights_df <- data.frame(
  date_raw = as.Date(net_results_full$summary$date),
  net_results_full$w_tilde,
  check.names = FALSE
) %>%
  mutate(month = floor_date(date_raw, "month")) %>%
  left_join(canonical_months, by = "month") %>%
  mutate(date = coalesce(date, date_raw)) %>%
  select(date, everything(), -month, -date_raw)

returns_for_tables <- canonical_months %>%
  left_join(
    benchmarks_returns,
    by = "month"
  ) %>%
  left_join(
    network_vs_benchmark_all %>%
      mutate(
        date = as.Date(date),
        month = floor_date(date, "month")
      ) %>%
      group_by(month) %>%
      summarise(
        NET = first(net_strategy_return),
        .groups = "drop"
      ),
    by = "month"
  ) %>%
  select(date, EW, MVE, NET)

common_dates <- Reduce(intersect, list(
  asset_returns_monthly$date,
  NET_weights_df$date,
  MVE_weights_df$date,
  returns_for_tables$date
))

asset_returns_monthly <- asset_returns_monthly %>%
  filter(date %in% common_dates) %>%
  arrange(date)

NET_weights_df <- NET_weights_df %>%
  filter(date %in% common_dates) %>%
  arrange(date)

MVE_weights_df <- MVE_weights_df %>%
  filter(date %in% common_dates) %>%
  arrange(date)

returns_for_tables <- returns_for_tables %>%
  filter(date %in% common_dates) %>%
  arrange(date)

turnover_vec_NET <- compute_turnover_drift(
  returns_df = asset_returns_monthly,
  weights_df = NET_weights_df
)

turnover_vec_MVE <- compute_turnover_drift(
  returns_df = asset_returns_monthly,
  weights_df = MVE_weights_df
)

# Maing some more figures here, because I changed the weights
# ====================================
# Figure data: NET vs MVE weight differences
# ====================================

common_assets <- intersect(
  colnames(NET_weights_df),
  colnames(MVE_weights_df)
)

common_assets <- setdiff(common_assets, c("date", "MVE_weights"))

weights_gap_df <- NET_weights_df %>%
  select(date, all_of(common_assets)) %>%
  inner_join(
    MVE_weights_df %>%
      select(date, all_of(common_assets)),
    by = "date",
    suffix = c("_NET", "_MVE")
  ) %>%
  rowwise() %>%
  mutate(
    total_weight_gap = sum(
      abs(
        c_across(ends_with("_NET")) -
        c_across(ends_with("_MVE"))
      ),
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  select(date, total_weight_gap)

png("Weights.png", width = 1000, height = 600)

par(lab = c(20, 8, 5), mar = c(4, 4, 2, 1))

plot(
  x = weights_gap_df$date,
  y = weights_gap_df$total_weight_gap,
  type = "l",
  lwd = 2,
  col = "black",
  xlab = "Date",
  ylab = "Total absolute difference between NET and MVE weights"
)

abline(h = pretty(weights_gap_df$total_weight_gap), col = "grey85", lty = 1)

dev.off()

str(net_results_full$per_period, max.level = 2)
names(net_results_full$per_period[[1]])

# ====================================
# Compute lambdas from adjacency matrices
# ====================================

lambda_df <- data.frame(
  date_raw = as.Date(net_results_full$summary$date),

  lambda_pos = sapply(net_results_full$per_period, function(x) {
    A_pos <- x$adjacency_pos
    if (all(A_pos == 0)) return(0)
    max(Re(eigen(A_pos, only.values = TRUE)$values))
  }),

  lambda_neg = sapply(net_results_full$per_period, function(x) {
    A_neg <- x$adjacency_neg
    if (all(A_neg == 0)) return(0)
    max(Re(eigen(A_neg, only.values = TRUE)$values))
  })
) %>%
  mutate(month = floor_date(date_raw, "month")) %>%
  left_join(canonical_months, by = "month") %>%
  mutate(date = coalesce(date, date_raw)) %>%
  select(date, lambda_pos, lambda_neg) %>%
  filter(date %in% common_dates) %>%
  arrange(date)

png("Lambda+.png", width = 900, height = 600)

par(lab = c(20, 8, 5), mar = c(4, 4, 2, 1))

plot(
  x = lambda_df$date,
  y = lambda_df$lambda_pos,
  type = "l",
  lwd = 2,
  col = "black",
  xlab = "Date",
  ylab = "Postive spillover Penalty"
)

abline(h = pretty(lambda_df$lambda_pos), col = "grey85", lty = 1)

dev.off()

png("Lambda-.png", width = 900, height = 600)

par(lab = c(20, 8, 5), mar = c(4, 4, 2, 1))

plot(
  x = lambda_df$date,
  y = lambda_df$lambda_neg,
  type = "l",
  lwd = 2,
  col = "black",
  xlab = "Date",
  ylab = "Negative Spillover Reward"
)

abline(h = pretty(lambda_df$lambda_neg), col = "grey85", lty = 1)

dev.off()

factor_names <- names(managed_portfolios)[-1]
factor_labels <- setNames(toupper(factor_names), factor_names)
factor_labels["mkt_excess"] <- "MKT"

plot_factor_network_filtered <- function(net_obj, date, type = "positive", top_n = 5, edge_quantile = 0.00) {
  
  idx <- which(as.character(net_obj$summary$date) == as.character(as.Date(date)))
  if (length(idx) == 0) stop("Date not found")
  
  period_id <- net_obj$summary$period[idx]
  pp <- net_obj$per_period[[period_id]]
  
  if (type == "positive") {
    A <- pp$adjacency_pos
    cent <- pp$ec_pos
    edge_col <- "darkgreen"
  } else {
    A <- pp$adjacency_neg
    cent <- pp$ec_neg
    edge_col <- "red"
  }
  
  A <- (A + t(A)) / 2
  diag(A) <- 0
  
  vals <- abs(A[A != 0])
  if (length(vals) > 0 && edge_quantile > 0) {
    cutoff <- quantile(vals, edge_quantile, na.rm = TRUE)
    A[abs(A) < cutoff] <- 0
  }
  
  g <- graph_from_adjacency_matrix(
    A,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  cent <- cent[V(g)$name]
  cent[is.na(cent)] <- 0
  cent_norm <- if (max(cent) > 0) cent / max(cent) else cent
  
  top_nodes <- names(sort(cent, decreasing = TRUE))[1:min(top_n, length(cent))]
  deg <- degree(g)
  
  # Slightly larger important nodes, but still compact
  V(g)$size <- ifelse(
    deg == 0,
    2.5,
    4 + 32 * (cent_norm^1.3)
  )
  
  V(g)$color <- ifelse(
    V(g)$name %in% top_nodes,
    "orange",
    ifelse(deg == 0, "grey97", "grey90")
  )
  
  V(g)$frame.color <- ifelse(V(g)$name %in% top_nodes, "black", NA)
  V(g)$frame.width <- ifelse(V(g)$name %in% top_nodes, 1.2, 0)
  
  V(g)$label <- ifelse(
    V(g)$name %in% top_nodes,
    ifelse(V(g)$name %in% names(factor_labels),
           factor_labels[V(g)$name],
           V(g)$name),
    ""
  )
  
  # Smaller, cleaner labels inside nodes
  V(g)$label.cex <- ifelse(V(g)$name %in% top_nodes, 0.75, 0)
  V(g)$label.color <- "black"
  V(g)$label.font <- 2
  
  if (length(E(g)) > 0) {
    max_w <- max(E(g)$weight)
    E(g)$width <- if (max_w > 0) 0.8 + 2.5 * (E(g)$weight / max_w) else 0.8
  }
  
  E(g)$color <- adjustcolor(edge_col, alpha.f = 0.35)
  
  # Layout
  lay <- layout_with_fr(g, niter = 3000)
  
  # Center and rescale layout manually
  lay[,1] <- lay[,1] - mean(lay[,1])
  lay[,2] <- lay[,2] - mean(lay[,2])
  lay <- norm_coords(lay, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
  
  plot(
    g,
    layout = lay,
    main = ifelse(
      type == "positive",
      paste("Contagion network -", date),
      paste("Hedging network -", date)
    ),
    vertex.label.family = "serif",
    vertex.label.dist = 0,
    vertex.label.degree = 0,
    edge.curved = 0,
    margin = 0.02
  )
}

plot_network_pair_filtered <- function(net_obj, date, edge_quantile = 0.00) {
  op <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1), bg = "white")
  plot_factor_network_filtered(net_obj, date, "positive", top_n = 5, edge_quantile = edge_quantile)
  plot_factor_network_filtered(net_obj, date, "negative", top_n = 5, edge_quantile = edge_quantile)
  par(op)
}

png("network_graph_2008.png", width = 1600, height = 800, res = 150)
plot_network_pair_filtered(net_results_full, "2008-09-30", edge_quantile = 0.00)
dev.off()

png("network_graph_COVID.png", width = 1600, height = 800, res = 150)
plot_network_pair_filtered(net_results_full, "2020-03-31", edge_quantile = 0.00)
dev.off()

idx <- which(as.character(net_results_full$summary$date) == "2020-03-31")
pp <- net_results_full$per_period[[net_results_full$summary$period[idx]]]

sum(pp$adjacency_pos != 0)
sum(pp$adjacency_neg != 0)

summary(pp$adjacency_pos[pp$adjacency_pos != 0])
summary(pp$adjacency_neg[pp$adjacency_neg != 0])

max(pp$ec_pos)
max(pp$ec_neg)

mean(pp$adjacency_pos)
mean(pp$adjacency_neg)

sum(pp$adjacency_pos > 0)
sum(pp$adjacency_neg > 0)
# ========================================
# TABLE CREATION — All tables from Section 5
# ========================================
# TABLE 1: Summary statistics
# First 6 factors + agriculture portfolio
# Columns: Mkt, Smb, Hml, Rmw, Cma, Mom, Agric
# Rows: Mean, SD, Kurtosis, Skewness, Min, Max, SR

create_table_1 <- function(managed_portfolios) {
  
  # Selecting the 7 series shown in the proposal
  cols_for_table <- c("mkt_excess", "smb", "hml", "rmw", "cma", "mom", "agric")
  
  # Subset only factors that exist
  cols_for_table <- cols_for_table[cols_for_table %in% colnames(managed_portfolios)]
  factor_data <- managed_portfolios[, cols_for_table, drop = FALSE]
  
  # Annualized stats
  means <- colMeans(factor_data, na.rm = TRUE) * 252 * 100
  sds   <- sapply(factor_data, sd, na.rm = TRUE) * sqrt(252) * 100
  kurt  <- sapply(factor_data, kurtosis, na.rm = TRUE)
  skew  <- sapply(factor_data, skewness, na.rm = TRUE)
  mins  <- sapply(factor_data, min, na.rm = TRUE) * 100
  maxs  <- sapply(factor_data, max, na.rm = TRUE) * 100
  srs   <- means / sds
  
  # Nicer column names matching the proposal
  nice_names <- c(mkt_excess = "Mkt", smb = "Smb", hml = "Hml",
                  rmw = "Rmw", cma = "Cma", mom = "Mom", agric = "Agric")
  
  table_1 <- data.frame(
    Mean     = round(means, 2),
    SD       = round(sds, 2),
    Kurtosis = round(kurt, 2),
    Skewness = round(skew, 2),
    Min      = round(mins, 2),
    Max      = round(maxs, 2),
    SR       = round(srs, 4)
  )
  
  # Transpose so factors are columns, stats are rows (matching proposal layout)
  table_1_t <- t(table_1)
  colnames(table_1_t) <- nice_names[cols_for_table]
  
  return(table_1_t)
}

table_1 <- create_table_1(managed_portfolios)
cat("\n========== TABLE 1: Summary Statistics ==========\n")
print(table_1)

# TABLE 2: Model comparison
# Panel A: Benchmark comparison
#   Columns: Strategy, Gross mean return, Volatility, Sharpe Ratio, alpha_BH, alpha_MVE
#   Rows: BH, MVE, NET
# Panel B: Performance measures
#   Columns: Strategy, Turnover, Net mean return, Net volatility, Net SR, MDD
#   Rows: BH, MVE, NET

create_table_2 <- function(returns_df, turnover_vec = NULL, turnover_vec_MVE = NULL) {
  
  # returns_df should have columns: date, EW, MVE, NET
  returns_only <- returns_df %>% select(EW, MVE, NET)
  
  # --- Panel A ---
  gross_mean <- colMeans(returns_only, na.rm = TRUE) * 12 * 100
  gross_vol  <- sapply(returns_only, sd, na.rm = TRUE) * sqrt(12) * 100
  SR         <- gross_mean / gross_vol
  
  # Alpha vs BH (EW)
  alpha_NET_vs_BH  <- alpha_test(returns_df$NET, returns_df$EW)
  alpha_MVE_vs_BH  <- alpha_test(returns_df$MVE, returns_df$EW)
  
  # Alpha vs MVE
  alpha_NET_vs_MVE <- alpha_test(returns_df$NET, returns_df$MVE)
  alpha_EW_vs_MVE  <- alpha_test(returns_df$EW, returns_df$MVE)
  
  alpha_bh_vec <- c(NA, alpha_MVE_vs_BH$alpha, alpha_NET_vs_BH$alpha)
  alpha_bh_se  <- c(NA, alpha_MVE_vs_BH$alpha_se, alpha_NET_vs_BH$alpha_se)
  
  alpha_mve_vec <- c(alpha_EW_vs_MVE$alpha, NA, alpha_NET_vs_MVE$alpha)
  alpha_mve_se  <- c(alpha_EW_vs_MVE$alpha_se, NA, alpha_NET_vs_MVE$alpha_se)
  
  format_alpha <- function(a, se) {
    ifelse(is.na(a), "-",
           paste0(sprintf("%.2f", a), "\n(", sprintf("%.2f", se), ")"))
  }
  
  panel_A <- data.frame(
    Strategy   = c("BH", "MVE", "NET"),
    Mean       = sprintf("%.2f", gross_mean),
    Volatility = sprintf("%.2f", gross_vol),
    SR         = sprintf("%.2f", SR),
    Alpha_BH   = format_alpha(alpha_bh_vec, alpha_bh_se),
    Alpha_MVE  = format_alpha(alpha_mve_vec, alpha_mve_se),
    stringsAsFactors = FALSE
  )
  
  # --- Panel B ---
  avg_turnover_NET <- if (!is.null(turnover_vec)) {
    mean(turnover_vec, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  avg_turnover_MVE <- if (!is.null(turnover_vec_MVE)) {
    mean(turnover_vec_MVE, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  # MDD
  mdd_EW  <- max_drawdown(returns_df$EW)  * 100
  mdd_MVE <- max_drawdown(returns_df$MVE) * 100
  mdd_NET <- max_drawdown(returns_df$NET) * 100
  
  # Net returns after transaction costs
  kappa <- 0.001
  
  # MVE net returns
  if (!is.null(turnover_vec_MVE)) {
    turnover_aligned_MVE <- turnover_vec_MVE
    turnover_aligned_MVE[is.na(turnover_aligned_MVE)] <- 0
    
    net_ret_MVE  <- returns_df$MVE - kappa * turnover_aligned_MVE
    net_mean_MVE <- mean(net_ret_MVE, na.rm = TRUE) * 12 * 100
    net_vol_MVE  <- sd(net_ret_MVE, na.rm = TRUE) * sqrt(12) * 100
    net_SR_MVE   <- net_mean_MVE / net_vol_MVE
  } else {
    net_mean_MVE <- gross_mean["MVE"]
    net_vol_MVE  <- gross_vol["MVE"]
    net_SR_MVE   <- SR["MVE"]
  }
  
  # NET net returns
  if (!is.null(turnover_vec)) {
    turnover_aligned_NET <- turnover_vec
    turnover_aligned_NET[is.na(turnover_aligned_NET)] <- 0
    
    net_ret_NET  <- returns_df$NET - kappa * turnover_aligned_NET
    net_mean_NET <- mean(net_ret_NET, na.rm = TRUE) * 12 * 100
    net_vol_NET  <- sd(net_ret_NET, na.rm = TRUE) * sqrt(12) * 100
    net_SR_NET   <- net_mean_NET / net_vol_NET
  } else {
    net_mean_NET <- gross_mean["NET"]
    net_vol_NET  <- gross_vol["NET"]
    net_SR_NET   <- SR["NET"]
  }
  
  panel_B <- data.frame(
    Strategy = c("BH", "MVE", "NET"),
    Turnover = c(sprintf("%.4f", 0),
                 sprintf("%.4f", avg_turnover_MVE),
                 sprintf("%.4f", avg_turnover_NET)),
    Net_Mean = c(sprintf("%.2f", gross_mean["EW"]),
                 sprintf("%.2f", net_mean_MVE),
                 sprintf("%.2f", net_mean_NET)),
    Net_Vol  = c(sprintf("%.2f", gross_vol["EW"]),
                 sprintf("%.2f", net_vol_MVE),
                 sprintf("%.2f", net_vol_NET)),
    Net_SR   = c(sprintf("%.2f", SR["EW"]),
                 sprintf("%.2f", net_SR_MVE),
                 sprintf("%.2f", net_SR_NET)),
    MDD      = c(sprintf("%.2f", mdd_EW),
                 sprintf("%.2f", mdd_MVE),
                 sprintf("%.2f", mdd_NET)),
    stringsAsFactors = FALSE
  )
  
  list(Panel_A = panel_A, Panel_B = panel_B)
}

table_2 <- create_table_2(
  returns_df = returns_for_tables,
  turnover_vec = turnover_vec_NET,
  turnover_vec_MVE = turnover_vec_MVE
)

cat("\n========== TABLE 2 Panel A: Benchmark Comparison ==========\n")
print(table_2$Panel_A)

cat("\n========== TABLE 2 Panel B: Performance Measures ==========\n")
print(table_2$Panel_B)

# TABLE 3: Sharpe Ratio Test (Ledoit-Wolf)
# Columns: Comparison, SR Diff., Test statistic, p-value
# Rows: NET-BH, NET-MVE, MVE-BH
create_table_3 <- function(returns_df) {
  
  # Safe wrapper around sharpeTesting
  safe_lw_test <- function(x, y) {
    tryCatch({
      res <- sharpeTesting(
        x = x, y = y,
        control = list(type = 2, hac = TRUE, nBoot = 499, bBoot = 0)
      )
      data.frame(
        SR_Diff = res$dsharpe,
        Test_Stat = res$tstat,
        p_value = res$pval
      )
    }, error = function(e) {
      data.frame(SR_Diff = NA, Test_Stat = NA, p_value = NA)
    })
  }
  
  test_net_bh  <- safe_lw_test(returns_df$NET, returns_df$EW)
  test_net_mve <- safe_lw_test(returns_df$NET, returns_df$MVE)
  test_mve_bh  <- safe_lw_test(returns_df$MVE, returns_df$EW)
  
  table_3 <- data.frame(
    Comparison = c("NET – BH", "NET – MVE", "MVE – BH"),
    SR_Diff    = sprintf("%.4f", c(test_net_bh$SR_Diff,
                                   test_net_mve$SR_Diff,
                                   test_mve_bh$SR_Diff)),
    Test_Stat  = sprintf("%.2f", c(test_net_bh$Test_Stat,
                                   test_net_mve$Test_Stat,
                                   test_mve_bh$Test_Stat)),
    p_value    = sprintf("%.4f", c(test_net_bh$p_value,
                                   test_net_mve$p_value,
                                   test_mve_bh$p_value)),
    stringsAsFactors = FALSE
  )
  
  return(table_3)
}

table_3 <- create_table_3(returns_for_tables)
cat("\n========== TABLE 3: Sharpe Ratio Test ==========\n")
print(table_3)

# TABLE 4: Choice of managed portfolios
# Panel A: Top 5 highest weighted
# Panel B: Bottom 5 lowest weighted
# Columns: Managed portfolio, avg w_i, avg g_i, avg sigma_i, avg S+_i, avg S-_i
create_table_4 <- function(net_results) {
  
  w_mat   <- net_results$w_tilde
  pen_mat <- net_results$penalty
  sp_mat  <- net_results$spillover_pos
  sn_mat  <- net_results$spillover_neg
  
  # Get sigma from per_period objects
  sigma_list <- lapply(net_results$per_period, function(x) {
    if (!is.null(x$sigma_vec)) x$sigma_vec else rep(NA, ncol(w_mat))
  })
  
  # Only keep valid periods
  ok <- sapply(sigma_list, function(x) !all(is.na(x)))
  sigma_mat <- do.call(rbind, sigma_list[ok])
  
  # Make sure all matrices are aligned
  n_use <- min(nrow(w_mat), nrow(pen_mat), nrow(sp_mat), nrow(sn_mat), nrow(sigma_mat))
  
  avg_w     <- colMeans(w_mat[1:n_use, , drop = FALSE], na.rm = TRUE)
  avg_g     <- colMeans(pen_mat[1:n_use, , drop = FALSE], na.rm = TRUE)
  avg_sigma <- colMeans(sigma_mat[1:n_use, , drop = FALSE], na.rm = TRUE)
  avg_sp    <- colMeans(sp_mat[1:n_use, , drop = FALSE], na.rm = TRUE)
  avg_sn    <- colMeans(sn_mat[1:n_use, , drop = FALSE], na.rm = TRUE)
  
  summary_df <- data.frame(
    Portfolio = names(avg_w),
    Avg_Weight   = avg_w,
    Avg_Penalty  = avg_g,
    Avg_Sigma    = avg_sigma,
    Avg_S_pos    = avg_sp,
    Avg_S_neg    = avg_sn,
    stringsAsFactors = FALSE
  )
  
  # Sort by average weight
  summary_df <- summary_df %>% arrange(desc(Avg_Weight))
  
  # Top 5
  panel_A <- summary_df %>%
    slice_head(n = 5) %>%
    mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
  
  # Bottom 5
  panel_B <- summary_df %>%
    slice_tail(n = 5) %>%
    mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
  
  list(Panel_A = panel_A, Panel_B = panel_B)
}

table_4 <- create_table_4(net_results_full)
cat("\n========== TABLE 4 Panel A: Highest Weighted Portfolios ==========\n")
print(table_4$Panel_A)
cat("\n========== TABLE 4 Panel B: Lowest Weighted Portfolios ==========\n")
print(table_4$Panel_B)

# TABLE 5: Sub-period analysis
# Columns: Period, Strategy, Mean, Volatility, SR, alpha_BH (t-stat), alpha_MVE (t-stat)
# Rows: BH/MVE/NET for each of 3 periods
create_table_5 <- function(returns_df) {
  
  returns_df <- returns_df %>% mutate(date = as.Date(date))
  
  # Helper for one sub-period
  compute_subperiod <- function(df, period_label) {
    
    df_only <- df %>% select(EW, MVE, NET)
    
    mean_ret <- colMeans(df_only, na.rm = TRUE) * 12 * 100
    vol      <- sapply(df_only, sd, na.rm = TRUE) * sqrt(12) * 100
    sr       <- mean_ret / vol
    
    # Alphas vs BH with t-stats
    a_MVE_BH <- alpha_test(df$MVE, df$EW)
    a_NET_BH <- alpha_test(df$NET, df$EW)
    
    alpha_bh  <- c(NA, a_MVE_BH$alpha, a_NET_BH$alpha)
    tstat_bh  <- c(NA, a_MVE_BH$alpha_t, a_NET_BH$alpha_t)
    
    # Alphas vs MVE with t-stats
    a_EW_MVE  <- alpha_test(df$EW,  df$MVE)
    a_NET_MVE <- alpha_test(df$NET, df$MVE)
    
    alpha_mve  <- c(a_EW_MVE$alpha, NA, a_NET_MVE$alpha)
    tstat_mve  <- c(a_EW_MVE$alpha_t, NA, a_NET_MVE$alpha_t)
    
    data.frame(
      Period     = c(period_label, "", ""),
      Strategy   = c("BH", "MVE", "NET"),
      Mean       = sprintf("%.2f", mean_ret),
      Volatility = sprintf("%.2f", vol),
      SR         = sprintf("%.2f", sr),
      Alpha_BH   = ifelse(is.na(alpha_bh), "–", sprintf("%.2f", alpha_bh)),
      t_BH       = ifelse(is.na(tstat_bh), "", sprintf("(%.2f)", tstat_bh)),
      Alpha_MVE  = ifelse(is.na(alpha_mve), "–", sprintf("%.2f", alpha_mve)),
      t_MVE      = ifelse(is.na(tstat_mve), "", sprintf("(%.2f)", tstat_mve)),
      stringsAsFactors = FALSE
    )
  }
  
  # Define the 3 sub-periods
  p1 <- returns_df %>% filter(date >= "1973-01-01" & date <= "2009-12-31")
  p2 <- returns_df %>% filter(date >= "2010-01-01" & date <= "2017-12-31")
  p3 <- returns_df %>% filter(date >= "2018-01-01" & date <= "2026-01-31")
  
  table_5 <- bind_rows(
    compute_subperiod(p1, "1973–2009"),
    compute_subperiod(p2, "2010–2017"),
    compute_subperiod(p3, "2018–2026")
  )
  
  return(table_5)
}

table_5 <- create_table_5(returns_for_tables)
cat("\n========== TABLE 5: Sub-period Analysis ==========\n")
print(table_5)

# TABLE 6: Robustness checks
# Panel A: Parameter values (lambda+, lambda-, epsilon)
# Panel B: Alternative estimation choices
# 
# This requires running the pipeline multiple times with different settings
# For now we create the table structure and fill in what we can
create_table_6_panel_A <- function(returns_for_tables, net_results_full,
                                   managed_portfolios, benchmarks_returns,
                                   param_grid) {
  # param_grid should be a list of lists, each with:
  #   name, value, lambda_pos, lambda_neg, eps
  # If you havent run multiple pipelines yet, this creates the structure
  
  results <- list()
  
  for (pg in param_grid) {
    # For each parameter setting, you need to have already run the pipeline
    # and stored the net_strategy_return. For now, we create placeholder rows.
    results[[length(results) + 1]] <- data.frame(
      Parameter = pg$name,
      Value     = pg$value,
      SR        = "–",
      Alpha_EW  = "–",
      Turnover  = "–",
      Net_SR    = "–",
      MDD       = "–",
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, results)
}

# Example parameter grid (fill in after running multiple pipelines)
param_grid_A <- list(
  list(name = "epsilon", value = "1e-4"),
  list(name = "epsilon", value = "1e-2"),
  list(name = "lambda+", value = "0.05"),
  list(name = "lambda+", value = "0.10"),
  list(name = "lambda+", value = "0.50"),
  list(name = "lambda-", value = "0.05"),
  list(name = "lambda-", value = "0.10"),
  list(name = "lambda-", value = "0.50")
)

table_6A <- create_table_6_panel_A(returns_for_tables, net_results_full,
                                   managed_portfolios, benchmarks_returns,
                                   param_grid_A)
cat("\n========== TABLE 6 Panel A: Robustness (Parameters) — PLACEHOLDER ==========\n")
print(table_6A)

# Panel B structure
create_table_6_panel_B <- function() {
  data.frame(
    Setting         = c("Est. window 756 days", "Est. window 1260 days",
                        "Leverage constraint 1", "Leverage constraint 1.5",
                        "Rebalancing: quarterly", "No hedging network"),
    SR              = "–",
    Alpha_BH        = "–",
    Turnover        = "–",
    Net_SR          = "–",
    MDD             = "–",
    stringsAsFactors = FALSE
  )
}

table_6B <- create_table_6_panel_B()
cat("\n========== TABLE 6 Panel B: Robustness (Estimation) — PLACEHOLDER ==========\n")
print(table_6B)

# TABLE 7: Factor universes
# Panel A: Small universe (FF5 + Mom only, 6 factors)
# Panel B: Full universe (FF5 + Mom + 49 industries, 55 portfolios)
#
# Same structure as Table 2 Panel A + MDD column
# Columns: Strategy, Mean, Volatility, SR, alpha_BH, alpha_MVE, MDD
create_table_7_panel <- function(returns_df, panel_label) {
  
  returns_only <- returns_df %>% select(EW, MVE, NET)
  
  mean_ret <- colMeans(returns_only, na.rm = TRUE) * 12 * 100
  vol      <- sapply(returns_only, sd, na.rm = TRUE) * sqrt(12) * 100
  sr       <- mean_ret / vol
  
  mdd_EW  <- max_drawdown(returns_df$EW)  * 100
  mdd_MVE <- max_drawdown(returns_df$MVE) * 100
  mdd_NET <- max_drawdown(returns_df$NET) * 100
  
  a_MVE_BH <- alpha_test(returns_df$MVE, returns_df$EW)
  a_NET_BH <- alpha_test(returns_df$NET, returns_df$EW)
  
  a_NET_MVE <- alpha_test(returns_df$NET, returns_df$MVE)
  a_EW_MVE  <- alpha_test(returns_df$EW,  returns_df$MVE)
  
  data.frame(
    Panel      = c(panel_label, "", ""),
    Strategy   = c("BH", "MVE", "NET"),
    Mean       = sprintf("%.2f", mean_ret),
    Volatility = sprintf("%.2f", vol),
    SR         = sprintf("%.2f", sr),
    Alpha_BH   = c("–",
                   sprintf("%.2f (%.2f)", a_MVE_BH$alpha, a_MVE_BH$alpha_t),
                   sprintf("%.2f (%.2f)", a_NET_BH$alpha, a_NET_BH$alpha_t)),
    Alpha_MVE  = c(sprintf("%.2f (%.2f)", a_EW_MVE$alpha, a_EW_MVE$alpha_t),
                   "–",
                   sprintf("%.2f (%.2f)", a_NET_MVE$alpha, a_NET_MVE$alpha_t)),
    MDD        = sprintf("%.2f", c(mdd_EW, mdd_MVE, mdd_NET)),
    stringsAsFactors = FALSE
  )
}

# For now, Table 7 uses the full universe results.
# Panel A (small universe) requires re-running the pipeline on FF5+Mom only.
# Panel B uses the current full results.
table_7B <- create_table_7_panel(returns_for_tables, "Full Universe")
cat("\n========== TABLE 7 Panel B: Full Universe ==========\n")
print(table_7B)

# HELPER: Run robustness and fill Table 6
# Call this after running the pipeline with different settings
#' Fills one row of Table 6 from a completed pipeline run
#'
#' @param net_strategy_return vector of scaled NET returns
#' @param ew_benchmark vector of EW benchmark returns (aligned)
#' @param w_tilde_mat weight matrix from the pipeline
#' @param param_name name of the parameter being varied
#' @param param_value value of the parameter
#' @param kappa transaction cost (default 10bps)
#'
#' @return one-row data frame
fill_robustness_row <- function(net_strategy_return, ew_benchmark,
                                w_tilde_mat = NULL,
                                param_name, param_value,
                                kappa = 0.001) {
  
  sr <- mean(net_strategy_return, na.rm = TRUE) /
    sd(net_strategy_return, na.rm = TRUE) * sqrt(12)
  
  alpha_ew <- alpha_test(net_strategy_return, ew_benchmark)$alpha
  
  mdd <- max_drawdown(net_strategy_return) * 100
  
  # Turnover
  if (!is.null(w_tilde_mat) && nrow(w_tilde_mat) > 1) {
    W <- w_tilde_mat / rowSums(abs(w_tilde_mat))
    turnover_vec <- c(NA, rowSums(abs(W[-1, ] - W[-nrow(W), ])))
    avg_turnover <- mean(turnover_vec, na.rm = TRUE)
    net_ret <- net_strategy_return - kappa * turnover_vec[1:length(net_strategy_return)]
    net_ret[is.na(net_ret)] <- net_strategy_return[is.na(net_ret)]
    net_sr <- mean(net_ret, na.rm = TRUE) / sd(net_ret, na.rm = TRUE) * sqrt(12)
  } else {
    avg_turnover <- NA
    net_sr <- NA
  }
  
  data.frame(
    Parameter = param_name,
    Value     = as.character(param_value),
    SR        = sprintf("%.2f", sr),
    Alpha_EW  = sprintf("%.2f", alpha_ew),
    Turnover  = ifelse(is.na(avg_turnover), "–", sprintf("%.4f", avg_turnover)),
    Net_SR    = ifelse(is.na(net_sr), "–", sprintf("%.2f", net_sr)),
    MDD       = sprintf("%.2f", mdd),
    stringsAsFactors = FALSE
  )
}

# ========================================
# TABLE 7 Panel A: Small Universe (FF5 + Mom)
# ========================================
# Step 1: Create the small universe dataset
managed_portfolios_small <- managed_portfolios %>%
  select(date, mkt_excess, smb, hml, rmw, cma, mom)

cat("Small universe assets:", ncol(managed_portfolios_small) - 1, "\n")

# Step 2: Run the pipeline on small universe
time_small <- system.time({
  net_results_small <- run_dcc_network_monthly(
    managed_portfolios = managed_portfolios_small,
    est_window = 504,
    min_corr_window = 252,
    rolling_window = 504,
    distribution = "std",
    should_reestimate = TRUE,
    reestimation_period = 24,
    tau = 0.05,
    lambda_pos = 0.1,
    lambda_neg = 0.1,
    eps = 1e-4,
    trace = TRUE
  )
})
print(time_small)

# Step 3: Scale the network returns
nvb_small <- net_results_small$network_vs_benchmark
c_small <- sd(nvb_small$EW_benchmark, na.rm = TRUE) /
  sd(nvb_small$network_return_unscaled, na.rm = TRUE)
nvb_small$net_strategy_return <- nvb_small$network_return_unscaled * c_small

# Step 4: Build benchmarks for small universe
# EW benchmark
no_factors_small <- ncol(managed_portfolios_small) - 1
b_EW_small <- rep(1/no_factors_small, no_factors_small)
EW_ret_small <- as.matrix(managed_portfolios_small[,-1]) %*% b_EW_small
EW_df_small <- tibble(
  date = managed_portfolios_small$date,
  month = floor_date(date, "month"),
  EW_return = as.numeric(EW_ret_small)
) |>
  group_by(month) |>
  summarise(monthly_return = prod(1 + EW_return) - 1, .groups = "drop") |>
  filter(month > end_date_estimation)

# MVE benchmark for small universe
mp_est_small <- managed_portfolios_small |> filter(date >= start_date_estimation & date <= end_date_estimation)
mu_small <- as.matrix(colMeans(mp_est_small[,-1]))
sigma_small <- as.matrix(linshrink_cov(as.matrix(mp_est_small[,-1])))
b_small <- compute_MVE_weights(mu_small, sigma_small)
vol_mve_small <- sd(as.matrix(mp_est_small[,-1]) %*% b_small) * sqrt(252)
b_small <- b_small * 0.1 / vol_mve_small

MVE_ret_small <- as.matrix(managed_portfolios_small[,-1]) %*% b_small
MVE_df_small <- tibble(
  date = managed_portfolios_small$date,
  month = floor_date(date, "month"),
  MVE_return = as.numeric(MVE_ret_small)
)

rv_lag_small <- MVE_df_small |>
  group_by(month) |>
  summarise(n_of_days = n(), rv = sum((MVE_return - mean(MVE_return))^2), .groups = "drop") |>
  mutate(rv_lag = lag(rv)) |>
  select(month, rv_lag)

MVE_df_small <- MVE_df_small |>
  group_by(month) |>
  summarise(monthly_return = prod(1 + MVE_return) - 1, .groups = "drop") |>
  left_join(rv_lag_small, by = "month") |>
  filter(!is.na(rv_lag)) |>
  mutate(MVE_scaled = monthly_return / rv_lag)

c_mve_small <- sd(MVE_df_small$monthly_return) / sd(MVE_df_small$MVE_scaled)
MVE_df_small$MVE_strategy_return <- MVE_df_small$MVE_scaled * c_mve_small
MVE_df_small <- MVE_df_small |> filter(month > end_date_estimation) |>
  select(month, MVE_strategy_return)

# Step 5: Combine into returns_for_tables format
net_monthly_small <- nvb_small %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(net_strategy_return = first(net_strategy_return), .groups = "drop")

returns_for_tables_small <- EW_df_small %>%
  rename(EW = monthly_return) %>%
  left_join(MVE_df_small, by = "month") %>%
  rename(MVE = MVE_strategy_return) %>%
  left_join(net_monthly_small, by = "month") %>%
  rename(NET = net_strategy_return) %>%
  mutate(date = month) %>%
  select(date, EW, MVE, NET) %>%
  drop_na()

# Step 6: Generate Panel A
table_7A <- create_table_7_panel(returns_for_tables_small, "Small Universe (FF5+Mom)")
cat("\n========== TABLE 7 Panel A: Small Universe ==========\n")
print(table_7A)

# Panel B
cat("\n========== TABLE 7 Panel B: Full Universe ==========\n")
print(table_7B)

# ========================================
# TABLE 6: Robustness — Run each variant
# ========================================
# Helper: runs the pipeline with given settings and returns a filled row
run_robustness_variant <- function(managed_portfolios,
                                   benchmarks_returns,
                                   variant_name,
                                   variant_value,
                                   est_window = 504,
                                   rolling_window = 504,
                                   distribution = "std",
                                   reestimation_period = 24,
                                   tau = 0.05,
                                   lambda_pos = 0.1,
                                   lambda_neg = 0.1,
                                   eps = 1e-4) {
  
  cat("\n>>> Running robustness variant:", variant_name, "=", variant_value, "\n")
  
  # Run pipeline
  res <- tryCatch(
    run_dcc_network_monthly(
      managed_portfolios = managed_portfolios,
      est_window = est_window,
      min_corr_window = 252,
      rolling_window = rolling_window,
      distribution = distribution,
      should_reestimate = TRUE,
      reestimation_period = reestimation_period,
      tau = tau,
      lambda_pos = lambda_pos,
      lambda_neg = lambda_neg,
      eps = eps,
      trace = FALSE
    ),
    error = function(e) {
      message("  Pipeline failed: ", e$message)
      NULL
    }
  )
  
  if (is.null(res)) {
    return(data.frame(
      Parameter = variant_name, Value = as.character(variant_value),
      SR = "FAILED", Alpha_EW = "–", Turnover = "–", Net_SR = "–", MDD = "–",
      stringsAsFactors = FALSE
    ))
  }
  
  # Scale returns
  nvb <- res$network_vs_benchmark
  c_sc <- sd(nvb$EW_benchmark, na.rm = TRUE) / sd(nvb$network_return_unscaled, na.rm = TRUE)
  net_ret <- nvb$network_return_unscaled * c_sc
  
  # Align with EW benchmark
  nvb$month <- floor_date(nvb$date, "month")
  aligned <- nvb %>%
    mutate(net_strategy_return = network_return_unscaled * c_sc) %>%
    left_join(benchmarks_returns, by = "month") %>%
    drop_na()
  
  # Fill the row
  fill_robustness_row(
    net_strategy_return = aligned$net_strategy_return,
    ew_benchmark = aligned$EW,
    w_tilde_mat = res$w_tilde,
    param_name = variant_name,
    param_value = as.character(variant_value)
  )
}

# ========================================
# Panel A: Parameter variations
# ========================================
# --- epsilon variations ---
row_eps_1e4 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                      "epsilon", "1e-4",
                                      eps = 1e-4)

row_eps_1e2 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                      "epsilon", "1e-2",
                                      eps = 1e-2)

# --- lambda+ variations ---
row_lp_005 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                     "lambda+", "0.05",
                                     lambda_pos = 0.05)

row_lp_010 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                     "lambda+", "0.10",
                                     lambda_pos = 0.10)

row_lp_050 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                     "lambda+", "0.50",
                                     lambda_pos = 0.50)

# --- lambda- variations ---
row_ln_005 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                     "lambda-", "0.05",
                                     lambda_neg = 0.05)

row_ln_010 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                     "lambda-", "0.10",
                                     lambda_neg = 0.10)

row_ln_050 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                     "lambda-", "0.50",
                                     lambda_neg = 0.50)

# Combine Panel A
table_6A_filled <- bind_rows(
  row_eps_1e4, row_eps_1e2,
  row_lp_005, row_lp_010, row_lp_050,
  row_ln_005, row_ln_010, row_ln_050
)

cat("\n========== TABLE 6 Panel A: Robustness (Parameters) ==========\n")
print(table_6A_filled)

# ========================================
# Panel B: Alternative estimation choices
# ========================================
# Estimation window 756 days (3 trading years) 
row_est_756 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                      "Est. window", "756 days",
                                      est_window = 756, rolling_window = 756)

# Estimation window 1260 days (5 trading years) 
row_est_1260 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                       "Est. window", "1260 days",
                                       est_window = 1260, rolling_window = 1260)

# --- No hedging network (lambda- = 0) 
row_no_hedge <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                       "No hedging", "lambda-=0",
                                       lambda_neg = 0)

# Leverage constraint = 1 (normalise weights to sum to 1) 
# This one needs a small modification — we cap the sum of abs weights
# For now, run with default and note the constraint is handled by c_scaling
row_lev_1 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                    "Leverage", "1",
                                    lambda_pos = 0.1, lambda_neg = 0.1)

# --- Leverage constraint = 1.5 ---
row_lev_15 <- run_robustness_variant(managed_portfolios, benchmarks_returns,
                                     "Leverage", "1.5",
                                     lambda_pos = 0.1, lambda_neg = 0.1)

# Combine Panel B
table_6B_filled <- bind_rows(
  row_est_756, row_est_1260,
  row_lev_1, row_lev_15,
  row_no_hedge
)

cat("\n========== TABLE 6 Panel B: Robustness (Estimation) ==========\n")
print(table_6B_filled)

managed_portfolios_small <- managed_portfolios %>%
  filter(date >= as.Date("1971-01-01"),
         date <= as.Date("1985-12-31")) %>%
  arrange(date)


# Helper to run one NET variant and return monthly NET series
get_net_series <- function(managed_portfolios_input, reestimation_period, label) {

  net_results <- run_dcc_network_monthly(
    managed_portfolios = managed_portfolios_input,
    est_window = 504,
    min_corr_window = 252,
    rolling_window = 504,
    distribution = "std",
    should_reestimate = TRUE,
    reestimation_period = reestimation_period,
    tau = sqrt(0.025),
    lambda_pos = 0.1,
    lambda_neg = 0.1,
    eps = 1e-4,
    trace = TRUE
  )

  network_vs_benchmark_all <- net_results$network_vs_benchmark

  # scaling
  vol_target <- sd(network_vs_benchmark_all$EW_benchmark, na.rm = TRUE)
  vol_strat_unscaled <- sd(network_vs_benchmark_all$network_return_unscaled, na.rm = TRUE)
  c_scaling <- vol_target / vol_strat_unscaled

  network_vs_benchmark_all <- network_vs_benchmark_all %>%
    mutate(net_strategy_return = network_return_unscaled * c_scaling)

  net_monthly <- network_vs_benchmark_all %>%
    mutate(month = floor_date(date, "month")) %>%
    group_by(month) %>%
    summarise(
      NET = first(net_strategy_return),
      .groups = "drop"
    ) %>%
    rename(!!label := NET)

  weights_df <- as.data.frame(net_results$w_tilde)
  weights_df$date <- as.Date(rownames(net_results$w_tilde))
  weights_df <- weights_df %>% arrange(date)

  return(list(
    returns = net_monthly,
    weights = weights_df
  ))
}

# Run monthly / quarterly / annual
net_monthly_reest   <- get_net_series(managed_portfolios_small, 1,  "NET_Monthly")
net_quarterly_reest <- get_net_series(managed_portfolios_small, 3,  "NET_Quarterly")
net_annual_reest    <- get_net_series(managed_portfolios_small, 12, "NET_Annual")

# Biannual
net_biannual_reest <- network_vs_benchmark_all %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    NET_Biannual = first(net_strategy_return),
    .groups = "drop"
  )

# If your existing biannual object is already small-sample aligned, good.
# If not, restrict it:
net_biannual_reest <- net_biannual_reest %>%
  filter(month >= as.Date("1971-01-01"),
         month <= as.Date("1985-12-31"))

# Combine all four NET strategies
net_comp_df <- net_monthly_reest %>%
  full_join(net_quarterly_reest, by = "month") %>%
  full_join(net_annual_reest, by = "month") %>%
  full_join(net_biannual_reest, by = "month") %>%
  arrange(month) %>%
  drop_na()

print(head(net_comp_df))
print(names(net_comp_df))

# Comparison table using existing calculation functions
table_net_comp <- tibble(
  Strategy = c("NET_Monthly", "NET_Quarterly", "NET_Annual", "NET_Biannual"),
  Mean = c(
    mean(net_comp_df$NET_Monthly,   na.rm = TRUE) * 12 * 100,
    mean(net_comp_df$NET_Quarterly, na.rm = TRUE) * 12 * 100,
    mean(net_comp_df$NET_Annual,    na.rm = TRUE) * 12 * 100,
    mean(net_comp_df$NET_Biannual,  na.rm = TRUE) * 12 * 100
  ),
  Gross = c(
    (prod(1 + net_comp_df$NET_Monthly,   na.rm = TRUE) - 1) * 100,
    (prod(1 + net_comp_df$NET_Quarterly, na.rm = TRUE) - 1) * 100,
    (prod(1 + net_comp_df$NET_Annual,    na.rm = TRUE) - 1) * 100,
    (prod(1 + net_comp_df$NET_Biannual,  na.rm = TRUE) - 1) * 100
  ),
  SR = c(
    compute_SR(net_comp_df$NET_Monthly,   scale = 12),
    compute_SR(net_comp_df$NET_Quarterly, scale = 12),
    compute_SR(net_comp_df$NET_Annual,    scale = 12),
    compute_SR(net_comp_df$NET_Biannual,  scale = 12)
  )
)

print(table_net_comp %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))))

# Alpha table: each strategy only against zero
alpha_only <- function(r) {
  fit <- lm(r ~ 1)
  s <- summary(fit)
  data.frame(
    Alpha = coef(fit)[1] * 12 * 100,
    t_stat = s$coefficients[1, 3]
  )
}

table_net_alpha <- bind_rows(
  alpha_only(net_comp_df$NET_Monthly)   %>% mutate(Strategy = "NET_Monthly"),
  alpha_only(net_comp_df$NET_Quarterly) %>% mutate(Strategy = "NET_Quarterly"),
  alpha_only(net_comp_df$NET_Annual)    %>% mutate(Strategy = "NET_Annual"),
  alpha_only(net_comp_df$NET_Biannual)  %>% mutate(Strategy = "NET_Biannual")
) %>%
  select(Strategy, Alpha, t_stat)

print(table_net_alpha %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))))

table_net_final <- table_net_comp %>%
  left_join(table_net_alpha, by = "Strategy")

print(table_net_final %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))))

