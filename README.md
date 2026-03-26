# network-analysis-volatility-management






# Object definitions:
# Data section (Kornel)

Important: 
- `start_date` - start date of the dataset (including)
- `end_date` - end date of the dataset (including)
- `factors_ff5_daily` - daily observations of Fama & French Five Factors returns
- `momentum_daily` - daily observations of Momentum factor returns 
- `industry_daily` - daily observations of industry managed portfolios returns
- `factors_joined_excess` - dataset containing: Fama French factors, momentum and industry managed portfolios after subtracting the risk free rate
- `managed_portfolios` - final dataset of factors and managed portfolios ready for the analysis 

For manipulation:
- `factors_joined` - joined Fama & French factors with momentum factor and industry managed portfolios (raw data)

### Benchmarks section (Kornel)

Important:
- `start_date_estimation` - start date of the estimation period
- `end_date_estimation` - end date of the estimation period
- `managed_portfolios_est` - subset of managed portfolio returns for estimation
- `mu` - sample mean
- `sigma` - Shrinkage covariance matrix (Shrinkage implementation from Ledoit and Wolf (2004) with identity matrix as the target matrix)
- `b` - in-sample optimal weights 
- `c` - scaling constant defined as sd(buy-and-hold)/sd(scaled MVE)
- `MVE_returns_df` - dataframe containing the returns of MVE strategy
- `no_of_factors` - number of factors for analysis
- `b_EW - weights` for equal weighted buy-and-hold
- `EW_returns_df` - dataframe containing returns of equal weighted buy-and-hold strategy
- `sharpe_ratios_benchmarks` - table with Sharpe Ratios for both benchmark strategies (annualized)
- `benchmarks_returns` - dataframe with returns of benchmarks for later performance comaprison

For manipulation:
- `sp_in_sample` - In-sample Sharpe Ratio
- `rv_lag_MVE` - data frame containing lagged realized variance
- `MVE_returns` - returns of MVE strategy before multiplying by c
- `EW_returns` - returns of equally weighted buy-and-hold strategy (without dates)

# DCC-NL Covariance Estimation
Important
- `managed_portfolios` – data frame containing daily returns of all factors and managed portfolios
- `dates_block` – vector of trading dates aligned across all series
- `init_window` – initial estimation window for univariate GARCH estimation, typically 504 trading days
- `min_corr_window` – minimum number of in-sample observations required for DCC estimation
- `Z_block` – matrix of standardised residuals from the univariate GARCH models
- `SIGMA_block` – matrix of conditional volatility forecasts from the univariate GARCH models
- `MU_block` – matrix of conditional mean forecasts from the univariate GARCH models
- `RESID_block` – matrix of residuals, defined as returns minus conditional mean forecasts
  
## Univariate GARCH Estimation
- `fit_garch_expanding_monthly()` – fits a univariate GARCH(1,1) model to one return series using an expanding window with monthly refits
- `run_all_garch_monthly()` – applies the monthly expanding-window GARCH estimation to all return series
- `mu` – one-step-ahead conditional mean forecast from the univariate GARCH model
- `sigma` – one-step-ahead conditional volatility forecast from the univariate GARCH model
- `sigma2` – one-step-ahead conditional variance forecast, equal to sigma^2
- `resid` – forecast error, defined as return minus conditional mean
- `z` – standardised residual, defined as resid / sigma
- `refit_date` – date on which the GARCH model is re-estimated

## GARCH Block Construction
- `build_garch_blocks()` – combines the list of univariate GARCH outputs into aligned matrices
- `Z_block` – T x N matrix of standardised residuals
- `SIGMA_block` – T x N matrix of conditional volatilities
- `MU_block` – T x N matrix of conditional means
- `RESID_block` – T x N matrix of residuals

## Monthly Forecasting Index
- `make_month_index()` – constructs the sequence of expanding in-sample windows and next-month forecast blocks
- `refit_month` – year-month identifier of the estimation month
- `refit_idx` – row index of the month-end estimation date
- `refit_date` – final trading day of the estimation month
- `insample_idx` – row indices used for in-sample DCC estimation
- `forecast_idx` – row indices corresponding to the next month forecast period
- `forecast_dates` – dates associated with the next month's forecast period
  
## Matrix Utilities
- `cov_to_cor()` – converts a covariance matrix into a correlation matrix
- `make_psd()` – forces a matrix to be positive semidefinite by truncating negative eigenvalues
- `safe_cor()` – computes a sample correlation matrix and applies PSD correction and normalisation
- `normalize_Q_to_R()` – converts the DCC intermediate matrix Q_t into a valid correlation matrix R_t

## Long-Run Correlation Target
- `S_target` – long-run correlation target used in the DCC recursion
- `safe_cor(Z_insample)` – benchmark correlation target based on the sample correlation of in-sample standardised residuals
- `build_nlshrink_target(Z_insample)` – nonlinear-shrinkage estimator of the covariance matrix, converted into a correlation     target for DCC
- `Sigma_nl` – nonlinear-shrunk covariance estimate based on the in-sample standardised residuals

The nonlinear-shrinkage target is used to stabilise the long-run correlation matrix in high-dimensional settings where the sample covariance or correlation matrix may be noisy or ill-conditioned.

## DCC Model Objects
- `a` – DCC news parameter, controlling the impact of the most recent shock
- `b` – DCC persistence parameter, controlling the influence of the previous DCC state
- `phi` – total persistence parameter, equal to a + b
- `S` – long-run correlation target used in the DCC recursion
- `Q_t` – intermediate covariance-like matrix in the DCC recursion
- `R_t` – conditional correlation matrix implied by Q_t
- `Q_prev` – lagged DCC state used in recursive updating
- `Q_T` – final in-sample DCC state used as the starting point for forecasting

## DCC Filtering and Estimation
- `dcc_filter()` – runs the DCC recursion over the standardised residual series and returns the sequence of Q_t and R_t
- `dcc_negloglik()` – computes the Gaussian negative log-likelihood for a given parameter vector (a, b)
- `estimate_dcc()` – estimates the DCC parameters (a, b) by minimising the negative log-likelihood subject to the standard DCC   constraints

## Parameter restrictions:
- a >= 0
- b >= 0
- a + b < 1

## Outputs from estimate_dcc():
- `a` – estimated DCC news parameter
- `b` – estimated DCC persistence parameter
- `S` – long-run correlation target used in estimation
- `Q_T` – final in-sample DCC state
- `opt` – optimisation output object
- `n_obs` – number of complete in-sample observations used in estimation

## Correlation Forecasting
- `forecast_dcc_correlations()` – produces multi-step-ahead forecasts of the DCC correlation matrices
- `dcc_fc` – list of future DCC states and correlation matrices
- `R_forecasts` – list of forecasted daily correlation matrices for the next month

## Covariance Reconstruction
- `build_cov_from_sigma_and_R()` – rebuilds daily covariance matrices using forecasted volatilities and forecasted correlations
- `D_t` – diagonal matrix of conditional volatilities at date t
- `H_t` – conditional covariance matrix at date t
- `H_forecasts` – list of daily covariance forecasts for the next month

## Monthly Covariance Aggregation
- `aggregate_monthly_cov()` – aggregates daily covariance forecasts into a monthly covariance matrix by summing valid daily      covariance forecasts
- `H_month` – monthly forecasted covariance matrix obtained by aggregating the daily covariance forecasts over the next month

## Master Function
- `run_dcc_monthly()` – executes the full monthly DCC forecasting procedure on expanding windows

For each monthly refit date, the function:
1. selects the in-sample block of standardised residuals
2. constructs the long-run correlation target
3. estimates the DCC parameters (a, b)
4. forecasts next-month daily correlation matrices
5. combines those correlation forecasts with the univariate GARCH volatility forecasts
6. aggregates the daily covariance forecasts into a monthly covariance matrix

Output Stored for Each Refit Month
- `refit_month` – month identifier of the estimation block
- `refit_date` – final trading day of the estimation month
- `insample_idx` – row indices used in DCC estimation
- `forecast_idx` – row indices used for forecasting
- `forecast_dates` – dates in the next-month forecast horizon
- `a` – estimated DCC news parameter
- `b` – estimated DCC persistence parameter
- `S` – long-run correlation target
- `Q_T` – final in-sample DCC state
- `R_forecasts` – forecasted daily correlation matrices
- `H_forecasts` – forecasted daily covariance matrices
- `H_month` – aggregated monthly covariance forecast
- `opt` – optimiser output from the DCC estimation step

Conceptual Flow
1. estimate univariate GARCH models for each return series
2. Extract standardised residuals and conditional volatility forecasts
3. build a stable long-run correlation target, optionally using nonlinear shrinkage
4. estimate DCC parameters on the in-sample standardised residuals
5. forecast next month's daily conditional correlations
6. Combine the forecasted correlations with the forecasted univariate volatilities
7. aggregate daily covariance forecasts into a monthly covariance matrix

### Network analysis (Zhang)
Important:
- `tau` - The threshold parameter (ranging from 0.1 to 0.5)
- `taus` – Sequence of threshold values used in the grid search (0.1 to 0.5)
- `grid_tau` - A list object storing the results of the grid search across 0.1 to 0.5
- `d` – Vector of inverse square roots of the diagonal elements of the precision matrix for scaling
- `adjacency_pos` - positive A+
- `adjacency_neg` - negative A-
- `ec_pos` - EC+
- `ec_neg` - EC-
- `spillover_pos` - S+
- `spillover_neg` - S-
- `partial_corr` - partial correlation
- `theta` - precision matrix
- `sigma_hat` – Forecasted covariance matrix
- `N` – Number of variables (dimension of the covariance/precision matrix)
- `eigenvalue_pos` – Eigenvalues of the positive adjacency matrix
- `eigenvector_pos` – Leading eigenvector of the positive adjacency matrix
- `ec_pos` – EC⁺, normalised to sum to 1.
- `eigenvalue_neg` – Eigenvalues of the negative adjacency matrix
- `eigenvector_neg` – Leading eigenvector of the negative adjacency matrix
- `ec_neg` – EC⁻, normalised to sum to 1.
- `sigma_vec` – Vector of individual asset volatilities (standard deviations)
- `temp_adjacency_pos` – Adjusted positive adjacency matrix with diagonal elements set to 0
- `temp_adjacency_neg` – Adjusted negative adjacency matrix with diagonal elements set to 0
- `init_spillover_pos` – Initial positive spillover effects before weighting by centrality
- `init_spillover_neg` – Initial negative spillover effects before weighting by centrality
- `spillover_pos` – S⁺
- `spillover_neg` – S⁻
- `adjacency_output` – Temporary object storing outputs of adjacency_matrix() during grid search
- `count_positive` – Number of positive links 
- `count_negative` – Number of negative links



# Libraries in the project
- `rlang`
- `tidyverse` - data manipulation after downloading
- `tidyfinance`
- `scales`
- `frenchdata` - downloading dataset from Kenneth Frenchs Library
- `dplyr`
- `moments` - library used to calculate higher order moments such as kurtosis and skewness
- `sandwich`
- `lmtest`
- `lubridate`
- `nlshrink` - used to estimate the covariance matrix with shrinkage from Ledoit & Wolf (2004)
- `rugarch`
- `igraph` - library used to make visualisations of the networks
