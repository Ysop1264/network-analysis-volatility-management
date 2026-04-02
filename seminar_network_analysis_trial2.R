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
# install.packages("xdcclarge")
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
    distribution = "norm",
    garch_order = c(1,1),
    arma_order = c(0,0)
){
  
  
  # Making sure of the form of returns and dates
  dates = as.Date(returns[[1]])
  returns <- as.numeric(returns[[2]])
  
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
    alpha_hat <- dcc_fit$a
    beta_hat  <- dcc_fit$b

    # Fallback if xdcclarge stores parameters differently
    if (is.null(alpha_hat) || length(alpha_hat) != 1 || !is.finite(alpha_hat) ||
        is.null(beta_hat)  || length(beta_hat)  != 1 || !is.finite(beta_hat)) {
    
      if (!is.null(dcc_fit$par) && length(dcc_fit$par) >= 2) {
        alpha_hat <- as.numeric(dcc_fit$par[1])
        beta_hat  <- as.numeric(dcc_fit$par[2])
      } else if (!is.null(dcc_fit$para) && length(dcc_fit$para) >= 2) {
        alpha_hat <- as.numeric(dcc_fit$para[1])
        beta_hat  <- as.numeric(dcc_fit$para[2])
      } else {
        message("xdcclarge returned object, but could not extract scalar alpha/beta.")
        message("Falling back to manual DCC estimation with NL target...")
        dcc_fit <- NULL
      }
    }

    if (!is.null(dcc_fit)) {
      if (!is.null(dcc_fit$C)) {
        target_used <- dcc_fit$C
      } else if (!is.null(target_C)) {
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
  # Indexing was off by one (added alternative below)
  # output_Q_list <- list()
  # for (i in 1:est_window){
    
  #   if(i == 1){
  #     output_Q_list[[i]] <- target_C
  #   } else {
  #     output_Q_list[[i]] <- (1-alpha-beta) * target_C + alpha * tcrossprod(as.matrix(residuals_dcc[i-1, ])) +
  #       beta * output_Q_list[[i-1]]
  #   }
    
  # }
  
  # Q <- output_Q_list[[est_window]]
  alpha <- as.numeric(alpha)
  beta  <- as.numeric(beta)
  target_C <- as.matrix(target_C)
  residuals_dcc <- as.matrix(residuals_dcc)

  if (length(alpha) != 1 || !is.finite(alpha)) stop("alpha must be one finite scalar")
  if (length(beta)  != 1 || !is.finite(beta)) stop("beta must be one finite scalar")
  if (!is.matrix(target_C) || nrow(target_C) != ncol(target_C)) stop("target_C must be square")
  if (ncol(residuals_dcc) != nrow(target_C)) stop("residuals_dcc columns must match target_C dimension")

  Q_prev <- target_C

  for(i in 1:est_window){
    zlag <- as.numeric(residuals_dcc[i, ])
    zlag[!is.finite(zlag)] <- 0
  
    if (length(zlag) != nrow(target_C)) {
    stop("zlag dimension mismatch in dcc_path")
    }
  
    zz <- tcrossprod(zlag)

    if (!all(dim(zz) == dim(target_C))) {
      stop("tcrossprod(zlag) dimension mismatch")
    }

    Q_prev <- (1 - alpha - beta) * target_C + alpha * zz + beta * Q_prev

    print(dim(Q_prev))
    print(class(Q_prev))
  
    Q_prev <- (Q_prev + t(Q_prev)) / 2  # keep symmetric
  }
  return(Q_prev)
}


# Function that aggregates the daily covariance forecasts to monthly covariance forecast
# Implements Eq. 6 in the proposal
# @param covariance_matricies list of daily forecasted matricies in the month
dcc_aggregate_daily_to_monthly <- function(covariance_matricies){
  monthly_covariance_matrix <- as.matrix(Reduce('+', covariance_matricies)/ length(covariance_matricies))
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
      print(length(dcc_est$alpha))
      print(dcc_est$beta)
      print(length(dcc_est$beta))
      print(class(dcc_est$target_C))
      print(dim(dcc_est$target_C))
      
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
      
      # Re-filter with frozen params but updated in-sample residuals for Q_T
      insample_df <- managed_portfolios[blk$insample_idx, ]
      
      garch_models <- tryCatch(
        estimate_all_univariate_garch(insample_df, est_window = est_window,
                                      distribution = distribution),
        error = function(e) {
          message("  GARCH refit failed, using frozen: ", e$message)
          NULL
        }
      )
      
      if (is.null(garch_models)) {
        garch_models <- frozen_garch_models
      }
      
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

# Test on smaller subset first
system.time({
 test_result <- run_dcc_network_monthly(
   managed_portfolios[1:1800, 1:8],
   est_window = 504,
   distribution = "norm",
   should_reestimate = FALSE,
   tau = 0.05,
   lambda_pos = 0.1,
   lambda_neg = 0.1,
   trace = TRUE
 )
})

test_result$summary
test_result$network_vs_benchmark
head(test_result$w_tilde)
head(test_result$penalty)

# Run for full sample
cat("Rows in full sample:", nrow(managed_portfolios), "\n")
cat("Number of assets/factors:", ncol(managed_portfolios) - 1, "\n")

time_full_pipeline <- system.time({
  net_results_full <- run_dcc_network_monthly(
    managed_portfolios,
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
# net_results_full$summary
# head(net_results_full$w_tilde)
# head(net_results_full$penalty)
# head(net_results_full$spillover_pos)
# head(net_results_full$spillover_neg)
# network_vs_benchmark_all


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
  rolling_SR_df$SR_NET[i-11] <- mean(combined_returns$net_strategy_return[(i-11):i]) / sd(combined_returns$net_strategy_return[(i-11):i]) * sqrt(12)
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
  
