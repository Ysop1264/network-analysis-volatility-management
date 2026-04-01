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
#install.packages("xdcclarge")
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





#====================================
# Settings for NET strategy
#=====================================


# Main settings
est_window = 504
#portfolios_set = c("mkt_excess", "smb", "hml", "rmw", "cma", "mom", "agric")
forecasting_period = 1000

# Restimation settings
should_reestimate_model = FALSE
reestimation_period = 252 # Now corresponds to a trading year

# GARCH(1,1) specification
garch_distribution_assumption = "norm"



# Setting smaller set for testing
#managed_portfolios <- managed_portfolios[1:2000, 1:6]
test_df <- managed_portfolios[1:1000, 1:6]





#'
#' Function that estimates GARCH(1,1) for each series
#' 
#' @param returns dataframe with dates as the first column, and managed portfolio returns as second
#' @param est_window estimation window length
#' @param distribution distribution assumption for GARCH(1,1) model
#' @param garch_order garch order
#' @param arma_model arma model for mean
#' @param scale_ret whether we multiply the returns by 100 for GARCH(1,1)
#' 
#' @return list with parameters, and estimated objects:mean, volatility, residuals, standardized residuals
#'
garch_estimate_single <- function(
    returns, 
    est_window = 504,
    distribution = "norm",
    garch_order = c(1,1),
    arma_order = c(0,0),
    scale_ret = FALSE
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
      include.mean = FALSE
    ),
    distribution.model = distribution
  )
  
  # Fitting the garch
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
    residuals_standard = c(residuals_st)#,
    #row.names = as.character(dates[1:est_window])
  )
  
  results_list <- list(
    parameters_model = parameters,
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
  lapply(colnames(returns_df)[-1], function(col_name) {
    
    # Create individual series data frame
    individual_series_df <- data.frame(
      date = as.Date(returns_df[[1]]),
      returns = returns_df[[col_name]]
    )
    
    garch_estimate_single(individual_series_df, est_window = est_window, distribution = distribution)
  })
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
  
  residuals_dcc <- do.call(cbind, 
                           lapply(univariate_garch_models, function(x) x$est_result$residuals_standard))
  colnames(residuals_dcc) <- colnames(returns_df)[-1]
  
  list <- list(
    H_initial = H,
    residuals_dcc_initial = residuals_dcc
  )
  
}

#'
#'Function that creates the target matrix with shrinakge method
#'
#'@param residuals_matrix matrix of standardized returns for DCC estimation
#'
#'@return C target matrix
create_target_matrix <- function(residuals_matrix){
  C <- cov2cor(nlshrink_cov(residuals_matrix))
}



#' 
#' Function that aggregates the daily covariance forecasts to monthly covariance forecast
#' 
#' @param covariance_matricies list of daily forecasted matricies in the month
#' 
#' 
dcc_aggregate_daily_to_monthly <- function(covariance_matricies){
  monthly_covariance_matrix <- as.matrix(reduce('+', covariance_matricies)/ length(covariance_matricies))
  
}

#' 
#' Function that estimates the path of the Q matrix in-sample period
#' 
#' @param alpha alpha from DCC estimate
#' @param beta beta from DCC estimate
#' @param target_C target unconditionall correlation matrix
#' @param residuals_dcc matrix with in-sample residuals
#' @param est_window length of estimation window
#' 
#' @return last estimated Q
#' 
dcc_path <- function(alpha, beta, target_C, residuals_dcc, est_window){
  
  output_Q_list <- list()
  for (i in 1:est_window){
    
    if(i == 1){
      output_Q_list[[i]] <- target_C
    } else {
      output_Q_list[[i]] <- (1-alpha-beta) * target_C + alpha * tcrossprod(as.matrix(residuals_dcc[i-1, ])) +
        beta * output_Q_list[[i-1]]
    }
    
    
    
  }
  
  Q <- output_Q_list[[est_window]]
}

#'
#' Function that forecasts the monthly covariance matrix from the DCC estimates
#'
#'@param alpha estimated alpha for DCC
#'@param beta estimated beta for DCC
#'@param final_Q last estimated Q from the in-sample period
#'@param target_C target unconditional correlation matrix from the in-sample
#'@param residuals_dcc matrix with in-sample residuals
#'@param no_days_month number of trading days in the considered month 
#'@param D_t matrix of forecasted volatilities from univariate GARCH
#'
#'@return estimated monthly covariance matrix
#'
dcc_forecast_monthly <- function(alpha, beta, final_Q, target_C, residuals_dcc, no_days_month, D_t){
  
  
  output_Q <- list()
  output_R <- list()
  output_H <- list()
  
  for(i in 1:no_days_month){
    
    if(i == 1){
      output_Q[[i]] <- (1 - alpha - beta) * target_C + 
        alpha * tcrossprod(as.matrix(residuals_dcc[nrow(residuals_dcc), ])) + beta * final_Q
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


#'
#'
#'Forecasts the univariate garch model for single return series
#'
#' @param returns dataframe with first column with dates and second column returns
#' @param est_window length of estimation window
#' @param fc_window forecasting window
#' @param garch_order order of GARCH model
#' @param arma_order order of ARMA model
#' @param distribution underlying distribution assumption for garch
#' @param scale_ret scaling returns or not (i will make better desc why)
#' 
#' @returns list with 1) parameters 2) forecasting results
#'
#'
forecast_single_univariate_garch <- function(
    returns, 
    est_window = 504,
    fc_window = 21,
    distribution = "norm",
    garch_order = c(1,1),
    arma_order = c(0,0),
    scale_ret = FALSE
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
      include.mean = FALSE
    ),
    distribution.model = distribution
  )
  
  # Fitting the garch
  fit <- ugarchfit(
    spec = spec_model, 
    data = returns, 
    out.sample = length(returns) - est_window, 
    solver = "hybrid"
    
  )
  
  fc <- ugarchforecast(
    fitORspec = fit,
    n.ahead = fc_window,
    n.roll = 0
    
  )
  
  #print("Length sigma" + length(sigma))
  fc_results <- data.frame(
    mu = c(fitted(fc)),
    sigma = c(sigma(fc))#,
    #row.names = as.character(dates[(est_window+1):length(returns)])
  )
  
  
  
  
  results_list <- list(
    fc_results = fc_results
  )
  
  return(results_list)
}



#'
#'Function that estimates the diagonal matrix of conditional volatilities based on univariate garch models
#'
#' @param returns dataframe with first column with dates and second column returns
#' @param est_window length of estimation window
#' @param fc_window forecasting window
#' @param distribution defning the underlying distributional assumption for garhch(1,1)
#' 
#' @return list of D_t matricies for forecasted period length 
#'
estimate_D_t <- function(returns_df,
                         est_window = 504, 
                         fc_window = 21, 
                         distribution = "norm"
){
  garch_results <- lapply(colnames(returns_df)[-1], function(col_name){
    # Create individual series data frame
    individual_series_df <- data.frame(
      date = as.Date(returns_df[[1]]),
      returns = returns_df[[col_name]]
    )
    
    result <- forecast_single_univariate_garch(individual_series_df, est_window = est_window, 
                                     fc_window = fc_window, distribution = distribution)
    return(result)
  })
  
  
  sigma_mat <- do.call(
    cbind,
    lapply(garch_results, function(x) x$fc_results$sigma)
  )
  
  colnames(sigma_mat) <- colnames(returns_df)[-1]
  
  D_list <- lapply(seq_len(nrow(sigma_mat)), function(t) {
    diag(sigma_mat[t, ])
  })
  
  return(list(
    garch_results = garch_results,
    sigma_mat = sigma_mat,
    D_list = D_list
  ))

}




