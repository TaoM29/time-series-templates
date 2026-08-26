# %% [markdown]
# # Daily Temperature Retrieval and Diagnostics
#
# This template retrieves daily mean air temperature from the Norwegian Meteorological Institute's Frost API, regularizes the observations, imputes short gaps, and demonstrates seasonal diagnostics. Generated values and plots are intentionally omitted; run the paired `Task2.R` script to reproduce them.
#
# Set the `FROST_CLIENT_ID` environment variable before requesting data. The example uses station `SN17850` from 2000 through 2023.
#
# %%

library(jsonlite)
library(tsibble)
library(dplyr)
library(tidyr)
library(lubridate)
library(imputeTS)
library(ggplot2)


# %% [markdown]
#
# ## Retrieve daily observations
#
# `get_temperature_tsibble()` requests the Frost element `mean(air_temperature P1D)`, keeps one observation per station and date, and returns a keyed daily `tsibble`.
#
# %%

get_temperature_tsibble <- function(client_id, station_id, start_date, end_date) {
  endpoint <- paste0("https://", client_id, "@frost.met.no/observations/v0.jsonld")
  elements <- "mean(air_temperature P1D)"

  url <- paste0(
    endpoint, "?",
    "sources=", station_id,
    "&referencetime=", start_date, "/", end_date,
    "&elements=", elements
  )
  
  response <- tryCatch(
    fromJSON(URLencode(url), flatten = TRUE),
    error = function(error) {
      stop("Unable to retrieve data from the Frost API: ", error$message, call. = FALSE)
    }
  )

  if (is.null(response$data)) {
    stop("The Frost API response did not contain observations.", call. = FALSE)
  }

  response$data %>%
    unnest(cols = observations) %>%
    transmute(
      Date = as.Date(referenceTime),
      sourceId,
      Temperature = value
    ) %>%
    distinct(Date, sourceId, .keep_all = TRUE) %>%
    as_tsibble(index = Date, key = sourceId)
}

client_id <- Sys.getenv("FROST_CLIENT_ID")
if (!nzchar(client_id)) {
  stop("Set FROST_CLIENT_ID before requesting data from frost.met.no.")
}

station_id <- "SN17850"
start_date <- "2000-01-01"
end_date <- "2023-12-31"

temperature_tsibble <- get_temperature_tsibble(client_id, station_id, start_date, end_date)


# %% [markdown]
#
# ## Regularize and impute the series
#
# Long early gaps can make a daily series unsuitable for seasonal analysis. The first helper keeps the suffix after the final gap longer than `max_gap`; the second fills missing calendar dates, records which values were absent, interpolates them, and removes leap days so every retained year has 365 observations.
#
# %%

print(has_gaps(temperature_tsibble))
print(count_gaps(temperature_tsibble))

find_acceptable_start <- function(tsibble_data, max_gap = 31) {
  gaps <- tsibble_data %>%
    arrange(Date) %>%
    group_by(sourceId) %>%
    mutate(gap_days = as.integer(Date - lag(Date))) %>%
    ungroup()
  
  acceptable_start_dates <- gaps %>%
    group_by(sourceId) %>%
    summarise(
      start_date = if (any(gap_days > max_gap, na.rm = TRUE)) {
        max(Date[gap_days > max_gap], na.rm = TRUE)
      } else {
        min(Date)
      },
      .groups = "drop"
    )

  tsibble_data %>%
    left_join(acceptable_start_dates, by = "sourceId") %>%
    filter(Date >= start_date) %>%
    select(-start_date) %>%
    as_tsibble(index = Date, key = sourceId)
}

process_temperature_tsibble <- function(tsibble_data, max_gap = 31) {
  limited_tsibble <- find_acceptable_start(tsibble_data, max_gap = max_gap)

  tsibble_regular <- limited_tsibble %>%
    group_by(sourceId) %>%
    complete(Date = seq(min(Date), max(Date), by = "day")) %>%
    ungroup() %>%
    as_tsibble(index = Date, key = sourceId)
  
  tsibble_imputed <- tsibble_regular %>%
    group_by(sourceId) %>%
    mutate(
      WasImputed = is.na(Temperature),
      Temperature = na_interpolation(Temperature)
    ) %>%
    ungroup()

  tsibble_imputed %>%
    filter(!(month(Date) == 2 & day(Date) == 29)) %>%
    as_tsibble(index = Date, key = sourceId)
}

processed_temperature_tsibble <- process_temperature_tsibble(temperature_tsibble)

density_data <- processed_temperature_tsibble %>%
  as_tibble() %>%
  mutate(ValueType = if_else(WasImputed, "Imputed", "Observed"))

ggplot(density_data, aes(x = Temperature, fill = ValueType)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Observed and Imputed Daily Temperatures",
    x = "Temperature",
    y = "Density",
    fill = "Value type"
  ) +
  theme_minimal()







# %% [markdown]
#
# ## Explore seasonality and dependence
#
# The processed values form a daily `ts` object with annual frequency 365. A long-lag ACF reveals recurring annual structure, the short-lag ACF focuses on persistence over four weeks, and the fixed-date scatter plot helps assess how the same point in the seasonal cycle changes across years.
#
# %%

temperature_ts <- ts(processed_temperature_tsibble$Temperature, 
                     frequency = 365)

ggplot(processed_temperature_tsibble, aes(x = Date, y = Temperature)) +
  geom_line() +
  labs(
    title = "Daily Mean Air Temperature",
    x = "Date",
    y = "Temperature"
  ) +
  theme_minimal()

# ACF plot up to 2000 lags (5.5 years)
acf_5_5_years <- acf(
  temperature_ts,
  lag.max = min(2000, length(temperature_ts) - 1),
  main = "ACF of Temperature (Up to 5.5 Years)"
)

# ACF plot up to 28 lags (4 weeks)
acf_4_weeks <- acf(
  temperature_ts,
  lag.max = min(28, length(temperature_ts) - 1),
  main = "ACF of Temperature (Up to 4 Weeks)"
)


# Scatter Plot for Specific Days Across Years (e.g., 1 October)
# Filter for October 1st entries
october_1st_data <- processed_temperature_tsibble %>%
  filter(format(Date, "%m-%d") == "10-01")


# Scatter plot of temperature on October 1st across years
ggplot(october_1st_data, aes(x = year(Date), y = Temperature)) +
  geom_point() +
  labs(
    title = "Temperature on 1 October Across Years",
    x = "Year",
    y = "Temperature"
  ) +
  theme_minimal()









