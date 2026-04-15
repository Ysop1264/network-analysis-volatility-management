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





# All factros and managed portfolios
asset_cols <- setdiff(names(managed_portfolios), "date")

# Only main factors
asset_cols_ff5_mom <- c("mkt_excess", "smb", "hml", "rmw", "cma", "mom")

# Industry managed portfolios
asset_cols_industries <- setdiff(names(industry_daily), "date")            


# set here whihc I want to do
returns_df <- managed_portfolios |>
  select(date, all_of(asset_cols))

# setting the lenght of the window and adding the observations (for main period there are non-missing tho)
window_days <- 252



# Getting the indicies function
get_month_end_indices <- function(dates) {
  split(seq_along(dates), format(dates, "%Y-%m")) |>
    purrr::map_int(max)
}

# setting the indicies for the figure 
month_end_idx <- get_month_end_indices(returns_df$date)
month_end_idx <- month_end_idx[month_end_idx >= window_days]

# computing the 'distribution'
corr_df <- purrr::map_dfr(month_end_idx, function(t_idx) {
  
  idx <- (t_idx - window_days + 1):t_idx
  R   <- as.matrix(returns_df[idx, asset_cols, drop = FALSE])
  
  # rolling correlation matrix
  C <- cor(R)
  
  # remove the diagonal entries of the correlation matrix
  diag(C) <- NA_real_
  
  avg_corr_i <- rowMeans(C, na.rm = TRUE)
  
  tibble(
    date     = returns_df$date[t_idx],
    q10      = as.numeric(quantile(avg_corr_i, 0.10, na.rm = TRUE)),
    q25      = as.numeric(quantile(avg_corr_i, 0.25, na.rm = TRUE)),
    q50      = as.numeric(quantile(avg_corr_i, 0.50, na.rm = TRUE)),
    q75      = as.numeric(quantile(avg_corr_i, 0.75, na.rm = TRUE)),
    q90      = as.numeric(quantile(avg_corr_i, 0.90, na.rm = TRUE)),
    mean     = mean(avg_corr_i, na.rm = TRUE)
  )
})








# constructing the plot
figure_1 <- ggplot(corr_df, aes(x = date)) +
  geom_ribbon(aes(ymin = q10, ymax = q90, fill = "80%"), alpha = 0.55) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = "50%"), alpha = 0.75) +
  geom_line(aes(y = q50, color = "Median"), linewidth = 1.0) +
  scale_fill_manual(
    values = c("50%" = "#5b9bd5", "80%" = "#c9d9ee"),
    name = NULL
  ) +
  scale_color_manual(
    values = c("Median" = "#1f4e79"),
    name = NULL
  ) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    x = "",
    y = "Average pairwise correlation",
    title = ""
  ) +
  guides(
    fill = guide_legend(order = 2),
    color = guide_legend(order = 1)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.minor = element_blank()
  )


# I am prining it for my R studio console (yall can comment it out)
print(figure_1)


# Saving this to the seminar folder 
ggsave(
  filename = "rolling_distribution_average_pairwise_correlation_wide.png",
  plot = figure_1,
  width = 14,
  height = 6,
  dpi = 300
)

