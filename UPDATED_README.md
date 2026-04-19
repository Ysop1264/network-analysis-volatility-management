# README — Network-Adjusted Volatility-Managed Portfolios Pipeline

This README documents the R script used to build benchmark portfolios, estimate a DCC-NL covariance process, construct network-adjusted portfolio weights, evaluate performance, and generate tables/figures.

---

## 1. Project purpose

The script builds and evaluates a **network-adjusted volatility-managed portfolio strategy** on Fama–French factors, momentum, and 49 industry portfolios.

The pipeline:

1. downloads daily factor and industry return data,
2. estimates univariate GARCH models for each series,
3. estimates dynamic correlations using a DCC framework with a nonlinear-shrinkage target,
4. aggregates daily covariance forecasts to the monthly level,
5. converts the covariance matrix into a signed network using partial correlations,
6. penalizes assets that are more exposed to positive spillovers and rewards those exposed to negative spillovers,
7. forms network-adjusted portfolio weights,
8. scales the resulting strategy to the volatility of an equal-weight benchmark,
9. compares the strategy to Buy-and-Hold (BH/EW) and MVE benchmarks.

---

## 2. Libraries used

- `rlang`
- `tidyverse`
- `tidyfinance`
- `scales`
- `frenchdata`
- `dplyr`
- `moments`
- `sandwich`
- `lmtest`
- `lubridate`
- `nlshrink`
- `rugarch`
- `igraph`
- `PeerPerformance`
- `xdcclarge`
- `tidyr`
- `purrr`

---

## 3. Data objects

### Main date controls
- `start_date`: start date of the downloaded sample.
- `end_date`: end date of the downloaded sample.

### Downloaded raw datasets
- `factors_ff5_daily`: daily Fama–French 5-factor data with market excess return and risk-free rate.
- `momentum_daily`: daily momentum factor returns.
- `industry_daily`: daily 49 industry portfolio returns.

### Combined datasets
- `factors_joined`: merged daily dataset containing FF5, momentum, and industry returns before excess-return conversion.
- `factors_joined_excess`: same merged dataset after subtracting the risk-free rate from all industry portfolios.
- `managed_portfolios`: final daily excess-return dataset used in the main analysis; excludes `risk_free`.

### Optional alternative universe
- `managed_portfolios_small`: appears twice in the script:
  - once as a **small asset universe** with only `mkt_excess`, `smb`, `hml`, `rmw`, `cma`, `mom`;
  - later overwritten as a **shorter time sample** from `1971-01-01` to `1985-12-31` for re-estimation frequency comparisons.

**Important:** the same object name is reused for two different purposes.

---

## 4. Core utility functions

### Portfolio utilities
- `compute_MVE_weights(mu, sigma)`: computes static mean-variance efficient weights as `solve(sigma, mu)` and normalizes them to sum to 1.
- `compute_SR(r)`: appears twice in the script.
  - first version: `mean(r) / sd(r)`;
  - later version: annualized Sharpe ratio using a `scale` argument.

### Plot helper
- `graph_annualized_returns(returns_df)`: base R plotting helper for annualized returns.

### Risk/performance helpers
- `annualized_mean(r, scale = 12)`: annualized mean return.
- `annualized_vol(r, scale = 12)`: annualized volatility.
- `max_drawdown(r)`: computes maximum drawdown from cumulative wealth.
- `compute_turnover(weights_df)`: simple turnover from changes in target weights between periods.
- `compute_net_returns(r, turnover, cost = 0.001)`: subtracts transaction cost times turnover.
- `performance_summary(r, turnover = NULL, scale = 12)`: summary of mean, volatility, Sharpe, drawdown, turnover, and net performance.
- `alpha_test(strategy_ret, benchmark_ret, lag = 3)`: OLS alpha regression with Newey–West standard errors.
- `sharpe_difference(r1, r2, scale = 12)`: difference in Sharpe ratios.
- `alpha_only(r)`: intercept-only regression used in the re-estimation comparison table.

---

## 5. Univariate GARCH block

### Single-series estimation
- `garch_estimate_single(...)`: fits a univariate GARCH(1,1) model to one return series.

Returns a list with:
- `parameters_model`: fitted GARCH parameters.
- `fit_object`: `ugarchfit` object.
- `est_result`: data frame containing:
  - `mu`: fitted conditional mean,
  - `sigma`: fitted conditional volatility,
  - `residuals`: raw residuals,
  - `residuals_standard`: standardized residuals.

**Implementation note:** returns are multiplied by 100 before fitting.

### Multi-series estimation
- `estimate_all_univariate_garch(returns_df, est_window, distribution)`: applies `garch_estimate_single()` to every managed portfolio series.

Returns a named list of fitted univariate GARCH models, one per asset/factor.

### DCC inputs built from GARCH outputs
- `creating_inputs_for_DCC(returns_df, univariate_garch_models)`: combines all fitted GARCH outputs into matrices.

Returns:
- `H_initial`: matrix of conditional variances (`sigma^2`) by date and asset.
- `residuals_raw`: matrix of raw residuals.
- `residuals_std`: matrix of standardized residuals.

---

## 6. Covariance target and matrix utilities

### Matrix stabilizers
- `make_psd(M, eps = 1e-8)`: forces a matrix to be positive semi-definite by replacing negative eigenvalues with `eps`.
- `cov_to_cor(Sigma, eps = 1e-10)`: converts a covariance matrix to a correlation matrix.

### Nonlinear-shrinkage target
- `create_target_matrix(residuals_matrix)`: estimates a nonlinear-shrunk covariance matrix from standardized residuals using `nlshrink_cov()`, converts it to correlation form, forces PSD, and sets diagonal elements to 1.

Returns:
- `C`: the long-run DCC correlation target.

---

## 7. DCC estimation and forecasting

### DCC parameter estimation
- `estimate_dcc_parameters(dcc_inputs, target_C = NULL)`: estimates DCC parameters `alpha` and `beta`.

Workflow:
1. tries `xdcclarge::cdcc_estimation(method = "NLS")`,
2. attempts to extract scalar `alpha` and `beta` from the returned object,
3. if that fails, falls back to a manual DCC likelihood optimization.

Returns:
- `alpha`: DCC news parameter.
- `beta`: DCC persistence parameter.
- `target_C`: correlation target used in estimation.
- `method`: either `"xdcclarge_NLS"` or `"manual_fallback"`.

### In-sample DCC recursion
- `dcc_path(alpha, beta, target_C, residuals_dcc, est_window)`: recursively updates `Q_t` through the in-sample period to recover the terminal state `Q_T`.

Returns:
- final `Q_prev`, i.e. the terminal in-sample DCC state.

### Aggregating daily covariances
- `dcc_aggregate_daily_to_monthly(covariance_matricies)`: sums daily forecast covariance matrices within a month.

Returns:
- `monthly_covariance_matrix`.

### Monthly covariance forecast
- `dcc_forecast_monthly(alpha, beta, final_Q, target_C, residuals_dcc, no_days_month, D_t)`: forecasts daily `Q_t`, transforms to `R_t`, reconstructs `H_t = D_t R_t D_t`, and aggregates daily matrices to monthly covariance.

Returns:
- monthly covariance matrix forecast `H_month`.

---

## 8. GARCH forecasting block

### One-step/multi-step forecast from fitted GARCH
- `forecast_from_fitted_garch(fit_object, fc_window = 21)`: produces `n.ahead = fc_window` forecasts from an already fitted `ugarchfit` model.

Returns a data frame with:
- `mu`: forecast conditional mean,
- `sigma`: forecast conditional volatility.

### Diagonal volatility matrices for DCC reconstruction
- `estimate_D_t(univariate_garch_models, fc_window = 21)`: forecasts all series and builds the diagonal volatility matrices used in covariance reconstruction.

Returns:
- `sigma_mat`: matrix of forecasted conditional volatilities.
- `D_list`: list of diagonal matrices `D_t` for each future trading day.

---

## 9. Network construction

### Adjacency matrices from partial correlations
- `adjacency_matrix(sigma_hat, tau = 0.05, ridge = 1e-6, use_correlation = TRUE)`: converts the monthly covariance matrix into a precision matrix, then into partial correlations.

Definitions used:
- `theta`: precision matrix, i.e. inverse of covariance or correlation matrix.
- `partial_corr`: matrix of pairwise partial correlations.
- `adjacency_pos`: positive network, keeps entries with partial correlation greater than `tau`.
- `adjacency_neg`: negative network, keeps entries with partial correlation less than `-tau`, stored in absolute value.

Returns:
- `partial_corr`
- `adjacency_pos`
- `adjacency_neg`

### Centrality measures
- `network_centrality(adjacency_pos, adjacency_neg)`: computes eigenvector centrality separately for the positive and negative adjacency matrices.

Returns:
- `ec_pos`: centrality in the contagion network.
- `ec_neg`: centrality in the hedging network.

### Spillover measures
- `spillovers(adjacency_pos, adjacency_neg, ec_pos, ec_neg, sigma_vec)`: computes spillover exposure for each asset.

Returns:
- `spillover_pos`: positive spillover exposure.
- `spillover_neg`: negative spillover exposure.

### Penalty and portfolio weights
- `network_penalty(spillover_pos, spillover_neg, lambda_pos = 0.1, lambda_neg = 0.1, eps = 1e-4)`: computes
  `g = max(eps, 1 + lambda_pos*S_pos - lambda_neg*S_neg)`.
- `network_adjusted_weights(...)`: computes raw network-adjusted weights
  `w_tilde = 1 / (sigma_vec * g)`.

Returns:
- `penalty`: vector `g`.
- `w_tilde`: raw network-adjusted weights.

### Portfolio return and scaling
- `network_portfolio_return(w_tilde, returns_next_month)`: computes the unscaled next-month portfolio return.
- `scaling_constant(f_net_unscaled, f_benchmark)`: scales the strategy to the volatility of a benchmark.
- `scaled_network_return(f_net_unscaled, c)`: applies the scaling constant.

---

## 10. Monthly indexing and re-estimation control

### Monthly forecast schedule
- `make_month_index(dates, min_obs = 252, rolling_window = NULL)`: builds the month-by-month estimation and forecast mapping.

Each list element contains:
- `refit_month`: estimation month identifier.
- `refit_idx`: row index of the refit date.
- `refit_date`: last trading day of the refit month.
- `insample_idx`: row indices used for estimation.
- `forecast_idx`: row indices for the next month.
- `forecast_dates`: trading dates of the forecast month.

### Re-estimation logic in the master function
The master function uses these internal objects for re-estimation control:
- `frozen_alpha`
- `frozen_beta`
- `frozen_target_C`
- `frozen_garch_models`
- `months_since_estimation`

If `should_reestimate = TRUE`, these objects are refreshed every `reestimation_period` months.

---

## 11. Master estimation function

### `run_dcc_network_monthly(...)`
This is the main pipeline.

For each eligible month it:
1. selects the in-sample window,
2. fits univariate GARCH models,
3. builds DCC inputs,
4. estimates the nonlinear-shrinkage correlation target,
5. estimates DCC parameters,
6. recovers the terminal in-sample `Q_T`,
7. forecasts next-month daily covariance matrices,
8. aggregates them into `H_month`,
9. constructs the positive/negative network,
10. computes spillovers and network-adjusted weights,
11. computes the realized next-month network portfolio return and equal-weight benchmark.

### Main arguments
- `managed_portfolios`: daily excess-return dataset.
- `est_window`: estimation window for GARCH fitting.
- `min_corr_window`: minimum sample size required to estimate DCC.
- `rolling_window`: rolling estimation length; if `NULL`, uses expanding window.
- `distribution`: innovation distribution for GARCH, e.g. `"norm"` or `"std"`.
- `should_reestimate`: whether to refresh the full GARCH+DCC system periodically.
- `reestimation_period`: number of months between re-estimations.
- `tau`: partial-correlation threshold for the adjacency matrix.
- `lambda_pos`, `lambda_neg`: network penalty parameters.
- `eps`: lower bound on the penalty.
- `trace`: whether to print progress and diagnostics.

### Output list
- `summary`: data frame containing valid forecast periods and dates.
- `w_tilde`: matrix of raw network-adjusted weights by month and asset.
- `penalty`: matrix of penalties by month and asset.
- `spillover_pos`: matrix of positive spillovers.
- `spillover_neg`: matrix of negative spillovers.
- `network_vs_benchmark`: data frame with:
  - `date`
  - `network_return_unscaled`
  - `EW_benchmark`
- `per_period`: full list of period-level objects.

### Important contents of each `per_period[[i]]`
- `refit_month`
- `refit_date`
- `insample_idx`
- `forecast_idx`
- `forecast_dates`
- `alpha`
- `beta`
- `H_month`
- `sigma_vec`
- `partial_corr`
- `adjacency_pos`
- `adjacency_neg`
- `ec_pos`
- `ec_neg`
- `spillover_pos`
- `spillover_neg`
- `penalty`
- `w_tilde`
- `realized_returns`
- `portfolio_return_unscaled`
- `EW_benchmark`

---

## 12. Main full-sample objects

### Network strategy
- `net_results_full`: main output from `run_dcc_network_monthly()` on the full universe.
- `network_vs_benchmark_all`: shorthand for `net_results_full$network_vs_benchmark`.
- `c_scaling`: scaling constant used to match the volatility of the network strategy to the EW benchmark.
- `net_strategy_return`: scaled network strategy return series.

### Benchmark estimation sample
- `start_date_estimation`: start of benchmark estimation period.
- `end_date_estimation`: end of benchmark estimation period.
- `managed_portfolios_est`: in-sample daily returns used to estimate MVE benchmark weights.

---

## 13. Benchmark construction

### MVE benchmark
- `mu`: in-sample mean vector of managed portfolios.
- `sigma`: shrinkage covariance matrix of managed portfolios from `linshrink_cov()`.
- `b`: static MVE portfolio weights, later scaled to target 10% annualized in-sample volatility.
- `vol_mve_in_sample`: in-sample annualized volatility of the unscaled MVE portfolio.
- `MVE_returns`: daily raw MVE portfolio return series.
- `MVE_returns_df`: monthly managed MVE return data frame.
- `rv_lag_MVE`: lagged monthly realized variance of the raw MVE portfolio.
- `MVE_scaled`: monthly MVE return divided by lagged realized variance.
- `c`: scaling constant used to match the volatility of the managed MVE strategy to the original monthly MVE return series.
- `MVE_strategy_return`: final scaled MVE monthly return.

### MVE weights through time
- `MVE_returns_weights_df`: monthly leverage scalar for the managed MVE strategy.
- `MVE_weights_mat`: monthly asset-level MVE weights, equal to the scalar monthly leverage multiplied by the static vector `b`.
- `MVE_weights_df`: data frame of monthly MVE asset weights aligned to canonical monthly dates.

### Equal-weight Buy-and-Hold benchmark
- `no_of_factors`: number of assets in the full universe.
- `b_EW`: static equal-weight vector.
- `EW_returns`: daily EW portfolio return series.
- `EW_returns_df`: monthly EW benchmark returns.

### Combined benchmark returns
- `benchmarks_returns`: data frame with columns:
  - `month`
  - `EW`
  - `MVE`

### Benchmark diagnostics
- `sharpe_ratios_benchmarks`: annualized Sharpe ratios for EW and MVE.

---

## 14. Aligned monthly objects used in evaluation

### Network/benchmark alignment
- `network_monthly`: monthly network strategy returns by month.
- `combined_returns`: merged monthly return data set with EW, MVE, and NET.
- `combined_returns_plot`: same as above but with complete cases only.
- `returns_for_tables`: final aligned performance-evaluation data frame with columns:
  - `date`
  - `EW`
  - `MVE`
  - `NET`

### Monthly asset returns and date alignment
- `asset_returns_monthly`: monthly compounded returns of every underlying asset/factor.
- `canonical_months`: helper data frame for aligning month identifiers with end-of-month dates.
- `common_dates`: common monthly dates across returns and weight objects.

### Network weights by month
- `NET_weights_df`: monthly raw network weights aligned to evaluation dates.

---

## 15. Turnover calculation

### Drift-based turnover
- `compute_turnover_drift(returns_df, weights_df, half_turnover = FALSE)`: computes turnover using drifting weights.

Mechanics:
1. start from previous target weights,
2. let them drift through realized asset returns,
3. compare the drifted weights to the new target weights,
4. sum absolute differences.

### Turnover vectors
- `turnover_vec_NET`: turnover series for the NET strategy.
- `turnover_vec_MVE`: turnover series for the MVE strategy.
- `EW_weights_df`: equal-weight holdings data frame used for turnover comparison.
- `turnover_vec_EW`: turnover series for the equal-weight benchmark.

---

## 16. Figures and figure objects

### Cumulative performance and rolling metrics
Files written by the script:
- `cumulative_wealth.jpeg`
- `rolling_sharpe.jpeg`
- `drawdown.jpeg`

Data objects:
- `cumulative_wealth_df`: cumulative wealth paths for EW, MVE, and NET.
- `rolling_SR_df`: rolling 12-month Sharpe ratios.
- `drawdown_df`: drawdown series.

### Network-adjustment mechanism figures
Files written by the script:
- `Weights.png`
- `Lambda+.png`
- `Lambda-.png`

Data objects:
- `weights_gap_df`: time series of total absolute differences between NET and MVE asset weights.
- `lambda_df`: time series of the largest eigenvalues of the positive and negative adjacency matrices.
  - `lambda_pos`: largest eigenvalue of the positive spillover network.
  - `lambda_neg`: largest eigenvalue of the negative spillover network.

### Network graph figures
Files written by the script:
- `network_graph_2008.png`
- `network_graph_COVID.png`
- `network_graph_1989.png` (opened but no `dev.off()` appears afterward in the pasted script)

Helper objects/functions:
- `factor_names`: raw asset/factor names.
- `factor_labels`: plotting labels.
- `plot_factor_network_filtered(...)`: plots a positive or negative network for a selected date.
- `plot_network_pair_filtered(...)`: plots the positive and negative network side by side.

---

## 17. Table construction functions

### Table 1
- `create_table_1(managed_portfolios)`: summary statistics for selected factors/portfolio series.
- `table_1`: output of Table 1.

### Table 2
- `create_table_2(returns_df, turnover_vec, turnover_vec_MVE, turnover_vec_EW)`: benchmark comparison and net-of-cost performance measures.
- `table_2`: output list with:
  - `Panel_A`: gross performance and alphas,
  - `Panel_B`: turnover, net returns, net Sharpe ratios, and maximum drawdown.

### Table 3
- `create_table_3(returns_df)`: Ledoit–Wolf Sharpe ratio difference tests.
- `table_3`: Sharpe-ratio comparison table.

### Table 4
- `create_table_4(net_results)`: average network weights, penalties, volatilities, and spillovers for the highest- and lowest-weighted portfolios.
- `table_4`: output list with `Panel_A` and `Panel_B`.

### Table 5
- `create_table_5(returns_df)`: sub-period performance analysis.
- `table_5`: output of Table 5.

### Table 6
- `create_table_6_panel_A(...)`: placeholder structure for parameter robustness.
- `create_table_6_panel_B()`: placeholder structure for estimation robustness.
- `fill_robustness_row(...)`: converts a completed robustness run into a table row.
- `run_robustness_variant(...)`: runs the main pipeline under alternative settings and formats the results.
- `table_6A`, `table_6B`: placeholder tables.
- `table_6A_filled`, `table_6B_filled`: filled robustness tables from variant runs.

### Table 7
- `create_table_7_panel(returns_df, panel_label)`: performance table for alternative asset universes.
- `table_7A`: small-universe results.
- `table_7B`: full-universe results.

---

## 18. Small-universe and robustness objects

### Small-universe run
- `net_results_small`: output from the NET pipeline using only FF5 + momentum.
- `nvb_small`: network-vs-benchmark returns for the small universe.
- `c_small`: scaling constant for small-universe NET.
- `EW_df_small`: equal-weight benchmark returns for the small universe.
- `mp_est_small`: estimation sample for small-universe MVE.
- `mu_small`, `sigma_small`, `b_small`: MVE inputs for the small universe.
- `MVE_df_small`: monthly managed MVE returns for the small universe.
- `returns_for_tables_small`: aligned performance data for small-universe Table 7.

### Robustness rows created in the script
Examples:
- `row_eps_1e4`, `row_eps_1e2`
- `row_lp_005`, `row_lp_010`, `row_lp_050`
- `row_ln_005`, `row_ln_010`, `row_ln_050`
- `row_est_756`, `row_est_1260`
- `row_no_hedge`
- `row_lev_1`, `row_lev_15`

---

## 19. Re-estimation frequency comparison objects

These objects compare monthly, quarterly, annual, and biannual re-estimation frequencies over a reduced sample.

### Helper function
- `get_net_series(managed_portfolios_input, reestimation_period, label)`: runs the NET pipeline, scales it to the EW benchmark, extracts monthly NET returns, and returns the raw results plus weight matrices.

Returns a list with:
- `returns`: monthly return series under the requested label.
- `weights`: monthly weight data frame.
- `raw`: original `run_dcc_network_monthly()` output.

### Re-estimation variants
- `net_monthly_reest`
- `net_quarterly_reest`
- `net_biannual_reest`
- `net_annual_reest`

### Comparison data frames
- `net_comp_df`: aligned monthly return series across re-estimation frequencies.
- `table_net_comp`: mean return, compounded gross return, and Sharpe ratio by frequency.
- `table_net_alpha`: alpha-only table by frequency.
- `table_net_final`: merged final comparison table.

---

## 20. Key implementation notes and caveats

1. **Returns are scaled by 100 inside the univariate GARCH estimation.**
   This affects fitted conditional means, volatilities, and residuals produced by `rugarch`, though the later pipeline works with the resulting objects consistently.

2. **`compute_SR()` is defined twice.**
   The second definition overwrites the first.

3. **`managed_portfolios_small` is reused for two different purposes.**
   This can be confusing when reading or extending the script.

4. **The DCC correlation target is estimated using nonlinear shrinkage on standardized residuals.**
   This is central to the DCC-NL implementation.

5. **The DCC estimation may come from `xdcclarge` or from a manual fallback optimizer.**
   The field `method` records which path was used.

6. **Re-estimation is genuinely periodic, not monthly by default.**
   When parameters are frozen, the previously estimated GARCH models, DCC parameters, and correlation target continue to be used until the next scheduled re-estimation.

7. **Turnover is currently computed without gross-exposure normalization.**
   The corresponding normalization lines are present but commented out.

8. **The script contains many diagnostic `print()` and `message()` calls.**
   These are useful for debugging but not part of the economic output objects.

9. **Some sections are placeholders or partial implementations.**
   In particular, parts of Table 6 are initially placeholder structures and are later filled using robustness reruns.
---

