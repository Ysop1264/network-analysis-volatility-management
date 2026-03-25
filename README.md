# network-analysis-volatility-management



## Object definitions:
### Data section (Kornel)

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

### Functions (Kornel & Yadhu)
- Functions are defined after data section and have corresponding documentation there

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

### Network analysis (Zhang)
Important:
- 'tau' - The threshold parameter (ranging from 0.1 to 0.5)
- 'taus' – Sequence of threshold values used in the grid search (0.1 to 0.5)
- 'grid_tau' - A list object storing the results of the grid search across 0.1 to 0.5
- 'd' – Vector of inverse square roots of the diagonal elements of the precision matrix for scaling
- 'adjacency_pos' - positive A+
- 'adjacency_neg' - negative A-
- 'ec_pos' - EC+
- 'ec_neg' - EC-
- 'spillover_pos' - S+
- 'spillover_neg' - S-
- 'partial_corr' - partial correlation
- 'theta' - precision matrix
- 'sigma_hat' – Forecasted covariance matrix
- 'N' – Number of variables (dimension of the covariance/precision matrix)
- 'eigenvalue_pos' – Eigenvalues of the positive adjacency matrix
- 'eigenvector_pos' – Leading eigenvector of the positive adjacency matrix
- 'ec_pos' – EC⁺, normalized to sum to 1.
- 'eigenvalue_neg' – Eigenvalues of the negative adjacency matrix
- 'eigenvector_neg' – Leading eigenvector of the negative adjacency matrix
- 'ec_neg' – EC⁻, normalized to sum to 1.
- 'sigma_vec' – Vector of individual asset volatilities (standard deviations)
- 'temp_adjacency_pos' – Adjusted positive adjacency matrix with diagonal elements set to 0
- 'temp_adjacency_neg' – Adjusted negative adjacency matrix with diagonal elements set to 0
- 'init_spillover_pos' – Initial positive spillover effects before weighting by centrality
- 'init_spillover_neg' – Initial negative spillover effects before weighting by centrality
- 'spillover_pos' – S⁺
- 'spillover_neg' – S⁻
- 'adjacency_output' – Temporary object storing outputs of adjacency_matrix() during grid search
- count_positive – Number of positive links 
- count_negative – Number of negative links
