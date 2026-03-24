# network-analysis-volatility-management



Object definitions:
Data section (Kornel)
Important: 
start_date - start date of the dataset (including)
end_date - end date of the dataset (including)
factors_ff5_daily - daily observations of Fama & French Five Factors returns
momentum_daily - daily observations of Momentum factor returns 
industry_daily - daily observations of industry managed portfolios returns
factors_joined_excess - dataset containing: Fama French factors, momentum and industry managed portfolios after subtracting the risk free rate
managed_portfolios - final dataset of factors and managed portfolios ready for the analysis 

For manipulation:
Factors_joined - joined Fama & French factors with momentum factor and industry managed portfolios (raw data)

Functions (Kornel & Yadhu)
Functions are defined after data section and have corresponding documentation there

Benchmarks section (Kornel)
Important:
start_date_estimation - start date of the estimation period
end_date_estimation - end date of the estimation period
managed_portfolios_est - subset of managed portfolio returns for estimation
mu - sample mean
sigma - Shrinkage covariance matrix (Shrinkage implementation from Ledoit and Wolf (2004) with identity matrix as the target matrix)
b - in-sample optimal weights 
c - scaling constant defined as sd(buy-and-hold)/sd(scaled MVE)
MVE_returns_df - dataframe containing the returns of MVE strategy
no_of_factors - number of factors for analysis
b_EW - weights for equal weighted buy-and-hold
EW_returns_df - dataframe containing returns of equal weighted buy-and-hold strategy
sharpe_ratios_benchmarks - table with Sharpe Ratios for both benchmark strategies (annualized)
benchmarks_returns - dataframe with returns of benchmarks for later performance comaprison

For manipulation:
sp_in_sample - In-sample Sharpe Ratio
rv_lag_MVE - data frame containing lagged realized variance
MVE_returns - returns of MVE strategy before multiplying by c
EW_returns - returns of equally weighted buy-and-hold strategy (without dates)
