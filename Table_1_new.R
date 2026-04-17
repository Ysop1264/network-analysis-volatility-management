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


# ============================================================
# Data download from th e tidy finance website (same as the main code) (factors and managed portoflios same as the first code)
# ============================================================
start_date <- as.Date("1971-01-01")
end_date   <- as.Date("2026-01-31")

factors_ff5_daily <- download_french_data("Fama/French 5 Factors (2x3) [Daily]")$subsets$data[[1]] |>
  mutate(
    date = ymd(date),
    across(c(RF, "Mkt-RF", SMB, HML, RMW, CMA), ~ as.numeric(.) / 100),
    .keep = "none"
  ) |>
  rename_with(str_to_lower) |>
  rename(mkt_excess = "mkt-rf", risk_free = rf) |>
  filter(date >= start_date & date <= end_date)

momentum_daily <- download_french_data("Momentum Factor (Mom) [Daily]")$subsets$data[[1]] |>
  mutate(
    date = ymd(date),
    across(Mom, ~ as.numeric(.) / 100),
    .keep = "none"
  ) |>
  rename_with(str_to_lower) |>
  filter(date >= start_date & date <= end_date)

industry_daily <- download_french_data("49 Industry Portfolios [Daily]")$subsets$data[[1]] |>
  mutate(
    date = ymd(date),
    across(-date, ~ as.numeric(.) / 100),
    .keep = "unused"
  ) |>
  rename_with(str_to_lower) |>
  filter(date >= start_date & date <= end_date)

factors_joined <- factors_ff5_daily |>
  left_join(momentum_daily, by = "date") |>
  left_join(industry_daily, by = "date")

factors_joined_excess <- factors_joined |>
  mutate(
    across(-c(date, mkt_excess, smb, hml, rmw, cma, risk_free, mom), ~ . - risk_free)
  )

managed_portfolios <- factors_joined_excess |>
  select(-risk_free)

returns_df <- managed_portfolios |>
  select(-date)

trading_days <- 252



# Defining functions for table computation (with annualization and changing to percentages)
ann_mean <- function(x, trading_days = 252) {
  mean(x) * trading_days * 100
}

ann_vol <- function(x, trading_days = 252) {
  sd(x) * sqrt(trading_days) * 100
}

ann_var <- function(x, trading_days = 252) {
  var(x) * trading_days * 100^2
}

ann_sharpe <- function(x, trading_days = 252) {
  s <- sd(x)
  mean(x) / s * sqrt(trading_days)
}




factor_stats <- tibble(
  factor      = names(returns_df),
  mean_return = sapply(returns_df, ann_mean, trading_days = trading_days),
  volatility  = sapply(returns_df, ann_vol, trading_days = trading_days),
  variance    = sapply(returns_df, ann_var, trading_days = trading_days),
  skewness    = sapply(returns_df, skewness),
  kurtosis    = sapply(returns_df, kurtosis),
  sharpe      = sapply(returns_df, ann_sharpe, trading_days = trading_days)
)



cs_summary <- function(x) {
  c(cs_mean = round(mean(x), 2),
    cs_sd   = round(sd(x), 2),
    p25     = round(as.numeric(quantile(x, 0.25, type = 7)), 2),
    median  = round(median(x), 2),
    p75     = round(as.numeric(quantile(x, 0.75, type = 7)), 2)
  )
}


panel_a <- bind_rows(
  tibble(
    Statistic = "Mean return",
    !!!as.list(cs_summary(factor_stats$mean_return))
  ),
  tibble(
    Statistic = "Sharpe ratio",
    !!!as.list(cs_summary(factor_stats$sharpe))
  )
)

panel_b <- bind_rows(
  tibble(
    Statistic = "Volatility",
    !!!as.list(cs_summary(factor_stats$volatility))
  ),
  tibble(
    Statistic = "Variance",
    !!!as.list(cs_summary(factor_stats$variance))
  ),
  tibble(
    Statistic = "Skewness",
    !!!as.list(cs_summary(factor_stats$skewness))
  ),
  tibble(
    Statistic = "Kurtosis",
    !!!as.list(cs_summary(factor_stats$kurtosis))
  )
)

