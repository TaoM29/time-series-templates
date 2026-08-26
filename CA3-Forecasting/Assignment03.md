# Electricity Demand Forecasting with Weather Covariates


This template combines daily private electricity consumption with temperature and global irradiation. It demonstrates time-zone-aware preprocessing, stationarity and dependence diagnostics, annual differencing, STL decomposition, rolling-origin ETS evaluation, Granger predictability tests, and ARIMAX holdout forecasting. Generated results and plots are intentionally omitted.

## Prepare electricity consumption

The source contains hourly readings with explicit UTC offsets. Converting them to `Europe/Oslo` before deriving calendar dates preserves the intended local day across daylight-saving transitions. The code identifies the longest contiguous suffix and aggregates consumption by day and consumer group.

```r

library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)

resolve_data_file <- function(filename, directory) {
  candidates <- c(
    filename,
    file.path(directory, filename),
    file.path("..", directory, filename)
  )
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop("Could not find ", filename, ". Run from the repository root or a template directory.")
  }

  existing[[1]]
}

# Read hourly consumption. read_csv2() respects decimal commas and the UTC
# offsets embedded in the ISO-8601 timestamps.
file_path <- resolve_data_file(
  "consumption_per_group_aas_hour.csv",
  "CA4-Classification"
)
consumption_hourly <- read_csv2(
  file_path,
  locale = locale(tz = "Europe/Oslo"),
  show_col_types = FALSE
)

# The course copy is already limited to Ås; this also makes a full Elhub
# municipality download safe to use under the same local filename.
if ("KOMMUNE" %in% names(consumption_hourly)) {
  consumption_hourly <- consumption_hourly %>%
    filter(KOMMUNE == "Ås")
}

# Select relevant columns
consumption_hourly <- consumption_hourly %>%
  select(STARTTID, FORBRUKSGRUPPE, VOLUM_KWH)


# Parse only if type inference did not already create a datetime column.
if (!inherits(consumption_hourly$STARTTID, "POSIXct")) {
  consumption_hourly <- consumption_hourly %>%
    mutate(STARTTID = parse_datetime(STARTTID, locale = locale(tz = "Europe/Oslo")))
}
consumption_hourly <- consumption_hourly %>%
  mutate(STARTTID = with_tz(STARTTID, "Europe/Oslo"))


# Extract unique dates from STARTTID and sort them
unique_dates <- sort(unique(as.Date(consumption_hourly$STARTTID)))

# Finding contiguous data start
date_diffs <- diff(unique_dates)

# Find where the difference is greater than 1 day, indicating a gap
gap_indices <- which(date_diffs > 1)
gaps <- data.frame(
  Start = unique_dates[gap_indices],
  End = unique_dates[gap_indices + 1] - 1
)

# Determine where the contiguous data starts
if (length(gap_indices) > 0) {
  contiguous_start <- unique_dates[gap_indices[length(gap_indices)] + 1]
} else {
  contiguous_start <- unique_dates[1]
}

# Reduce the data to daily sums for each type of consumer
daily_data <- consumption_hourly %>%
  group_by(Date = as.Date(STARTTID), FORBRUKSGRUPPE) %>%
  summarise(VOLUM_KWH = sum(VOLUM_KWH, na.rm = TRUE), .groups = "drop") %>%
  filter(Date >= contiguous_start) %>%
  ungroup()


print("Gaps in the dataset:")
print(gaps)
print(paste("Contiguous data starts from:", contiguous_start))

```

### Inspect daylight-saving transitions

In Oslo, the spring transition omits one local hour and the autumn transition repeats one. The following slices make those boundary days easy to inspect before daily aggregation.

```r

# Check for missing data in March 27, 2022
march_data <- consumption_hourly %>%
  filter(STARTTID >= as.POSIXct("2022-03-27 00:00:00", tz = "Europe/Oslo") & 
         STARTTID < as.POSIXct("2022-03-28 00:00:00", tz = "Europe/Oslo"))

# Check for duplicates in October 30, 2022
october_data <- consumption_hourly %>%
  filter(STARTTID >= as.POSIXct("2022-10-30 00:00:00", tz = "Europe/Oslo") & 
         STARTTID < as.POSIXct("2022-10-31 00:00:00", tz = "Europe/Oslo"))

# Print the relevant time sections
cat("Data around March 27, 2022 (expect missing 02:00):\n")
print(march_data)

cat("\nData around October 30, 2022 (expect duplicates at 02:00):\n")
print(october_data, n=13)

```



## Load weather covariates

The annual spreadsheets contain daily average air temperature (`LT`) and global irradiation (`GLOB`) from the NMBU meteorological station. `process_data_year()` standardizes column names and accepts either native spreadsheet dates or common textual date formats.

```r

library(readxl)

# Function to read and process each year's data
process_data_year <- function(year) {
  file_name <- resolve_data_file(
    paste0("Aas dogn ", year, ".xlsx"),
    "CA3-Forecasting"
  )
  
  # Read the Excel file
  data <- read_excel(file_name)

  # Extract relevant columns based on the actual names
  data_processed <- data %>%
    select(DATO, LT, GLOB) %>% 
    rename(
      date = DATO,
      avg_temp = LT,
      global_irradiation = GLOB
    )

  # readxl may return native dates or character values depending on the sheet.
  if (inherits(data_processed$date, c("Date", "POSIXt"))) {
    data_processed$date <- as.Date(data_processed$date)
  } else {
    data_processed$date <- as.Date(parse_date_time(
      as.character(data_processed$date),
      orders = c("ymd HMS", "dmy HMS", "mdy HMS", "ymd", "dmy", "mdy"),
      quiet = TRUE
    ))
  }
  
  return(data_processed)
}

# List of years to process
years <- 2017:2024 

# Combine data from all years
all_data <- bind_rows(lapply(years, process_data_year)) %>%
  arrange(date)

# View the combined data, to ensure the starting point (2017) and ending point (2024)
print(head(all_data))
print(tail(all_data))

```

### Verify leap-day handling

The retained dates are checked explicitly for 29 February in the leap years represented by the source files.

```r

# Check for leap years and February 29 entries
check_leap_years <- function(data) {
  leap_years <- c(2020, 2024)
  
  # Filter the dataset for leap years and look for February 29
  leap_day_data <- data %>%
    filter(year(date) %in% leap_years, month(date) == 2, day(date) == 29)
  
  return(leap_day_data)
}

# Call the function with your combined data
leap_day_entries <- check_leap_years(all_data)
print(leap_day_entries)

```

### Impute missing weather values

Linear interpolation uses neighboring dates for internal gaps, while `rule = 2` carries the nearest observed endpoint when a file begins or ends with missing values. This simple method is appropriate for short gaps but should be reconsidered for long missing intervals.

```r

library(zoo)  

# Check for missing values in avg_temp and global_irradiation
missing_data_summary <- all_data %>%
  summarise(
    missing_avg_temp = sum(is.na(avg_temp)),
    missing_global_irradiation = sum(is.na(global_irradiation))
  )

print(missing_data_summary)

# Impute missing values using linear interpolation for global irradiation
all_data <- all_data %>%
  mutate(
    avg_temp = na.approx(avg_temp, x = date, na.rm = FALSE, rule = 2),
    global_irradiation = na.approx(
      global_irradiation,
      x = date,
      na.rm = FALSE,
      rule = 2
    )
  )

# Confirm that there are no more missing values
missing_data_summary_after <- all_data %>%
  summarise(
    missing_avg_temp = sum(is.na(avg_temp)),
    missing_global_irradiation = sum(is.na(global_irradiation))
  )

print(missing_data_summary_after)



```

## Merge consumption and weather

The modeling table uses private consumption only, joins on the common daily range, and removes leap days so annual seasonal operations use a consistent 365-day cycle. Keeping one row per day avoids repeating weather values across hourly observations or mixing consumer groups in a single response series.

```r

# Keep one daily target series. Mixing consumer groups or repeating daily weather
# across hourly rows would invalidate the autocorrelation and forecast models.
private_daily <- daily_data %>%
  filter(FORBRUKSGRUPPE == "Privat")

start_date_longest <- max(min(private_daily$Date), min(all_data$date))
end_date_longest <- min(max(private_daily$Date), max(all_data$date))

merged_data_final <- private_daily %>%
  filter(Date >= start_date_longest, Date <= end_date_longest) %>%
  inner_join(
    all_data %>%
      filter(date >= start_date_longest, date <= end_date_longest),
    by = c("Date" = "date")
  ) %>%
  filter(!(month(Date) == 2 & day(Date) == 29)) %>%
  arrange(Date)

merged_data_clean <- merged_data_final %>%
  filter(if_all(c(VOLUM_KWH, avg_temp, global_irradiation), ~ !is.na(.x)))

print(paste("Merged daily observations:", nrow(merged_data_final)))
print(colSums(is.na(merged_data_clean)))


```


## Explore the aligned series

### Visualize the inputs

The hourly consumption plot is useful for checking consumer-group behavior and timestamp continuity. The daily weather plots expose long-term trends, annual seasonality, and remaining anomalies before modeling.

```r

# Hourly energy consumption by consumer group
ggplot(consumption_hourly, aes(x = STARTTID, y = VOLUM_KWH, color = FORBRUKSGRUPPE)) +
  geom_line() +
  labs(title = "Hourly Energy Consumption",
       x = "Date and Time",
       y = "Energy Consumption (kWh)",
       color = "Consumption Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Visualization for `all_data` - Daily Average Temperature
ggplot(all_data, aes(x = date, y = avg_temp)) +
  geom_line(color = "blue") +
  labs(title = "Daily Average Temperature",
       x = "Date",
       y = "Average Temperature (°C)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Visualization for `all_data` - Daily Global Irradiation
ggplot(all_data, aes(x = date, y = global_irradiation)) +
  geom_line(color = "orange") +
  labs(title = "Daily Global Irradiation",
       x = "Date",
       y = "Global Irradiation (W/m²)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


```

### Check stationarity

The ADF and KPSS tests use opposite null hypotheses: ADF tests for a unit root, while KPSS tests level stationarity. Reading them together helps distinguish clear stationarity from conflicting evidence caused by trend or seasonality. All tests below use the aligned daily modeling series.

```r

library(urca)
library(tseries)

# Private consumption
adf_test_volum <- ur.df(merged_data_clean$VOLUM_KWH, type = "drift", lags = 1)
summary(adf_test_volum)
kpss_test_volum <- kpss.test(merged_data_clean$VOLUM_KWH)
print(kpss_test_volum)

# Average temperature
adf_test_temp <- ur.df(merged_data_clean$avg_temp, type = "drift", lags = 1)
summary(adf_test_temp)
kpss_test_temp <- kpss.test(merged_data_clean$avg_temp)
print(kpss_test_temp)

# Global irradiation
adf_test_irrad <- ur.df(merged_data_clean$global_irradiation, type = "drift", lags = 1)
summary(adf_test_irrad)
kpss_test_irrad <- kpss.test(merged_data_clean$global_irradiation)
print(kpss_test_irrad)

```

### Inspect cross-correlation and serial dependence

The contemporaneous correlation matrix summarizes linear association among consumption, temperature, and irradiation. ACF/PACF plots then examine within-series dependence at short horizons and across roughly two years.

```r

library(forecast)

# Calculate cross-correlations
correlation_matrix <- merged_data_clean %>%
  select(VOLUM_KWH, avg_temp, global_irradiation) %>%
  cor(use = "complete.obs")

# Print the correlation matrix
print(correlation_matrix)




# Function to plot ACF and PACF
plot_acf_pacf <- function(data, variable, title_suffix, max_lag) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(1, 2))

  lag_limit <- min(max_lag, nrow(data) - 1)
  acf(data[[variable]], main = paste("ACF of", variable, title_suffix), lag.max = lag_limit)
  pacf(data[[variable]], main = paste("PACF of", variable, title_suffix), lag.max = lag_limit)
}

# Short-term analysis (up to 14 days)
plot_acf_pacf(merged_data_clean, "avg_temp", "(Short-term)", max_lag = 14)
plot_acf_pacf(merged_data_clean, "global_irradiation", "(Short-term)", max_lag = 14)
plot_acf_pacf(merged_data_clean, "VOLUM_KWH", "(Short-term)", max_lag = 14)

# Long-term analysis (minimum 2 years, assuming daily data ~ 730 days)
plot_acf_pacf(merged_data_clean, "avg_temp", "(Long-term)", max_lag = 730)
plot_acf_pacf(merged_data_clean, "global_irradiation", "(Long-term)", max_lag = 730)
plot_acf_pacf(merged_data_clean, "VOLUM_KWH", "(Long-term)", max_lag = 730)



```

### Interpreting ACF and PACF

Inspect the short-lag plots for immediate persistence and weekly structure, and the long-lag plots for slowly decaying or annually repeating dependence. A sharp PACF cutoff can suggest an autoregressive order, while a sharp ACF cutoff can suggest a moving-average order. Because outputs are stripped, run the paired script before recording series-specific conclusions.

## Remove annual seasonality

### Annual differencing

Subtracting each value from the observation 365 days earlier targets recurring yearly structure. Comparing the differenced ACF/PACF with the original plots shows which short-term dependencies remain.

```r

seasonal_lag <- 365

# Perform seasonal differencing for each variable
diff_avg_temp <- diff(merged_data_clean$avg_temp, lag = seasonal_lag)
diff_global_irradiation <- diff(merged_data_clean$global_irradiation, lag = seasonal_lag)
diff_VOLUM_KWH <- diff(merged_data_clean$VOLUM_KWH, lag = seasonal_lag)

  
# Function to plot ACF and PACF for short and long term
plot_acf_pacf_diff <- function(data, variable_name, short_lag, long_lag) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 2))

  short_lag <- min(short_lag, length(data) - 1)
  long_lag <- min(long_lag, length(data) - 1)
  
  # Short-term ACF and PACF
  acf(data, main = paste("ACF of", variable_name, "(Short-term, lag =", short_lag, ")"), lag.max = short_lag)
  pacf(data, main = paste("PACF of", variable_name, "(Short-term, lag =", short_lag, ")"), lag.max = short_lag)
  
  # Long-term ACF and PACF
  acf(data, main = paste("ACF of", variable_name, "(Long-term, lag =", long_lag, ")"), lag.max = long_lag)
  pacf(data, main = paste("PACF of", variable_name, "(Long-term, lag =", long_lag, ")"), lag.max = long_lag)
  
}


# Plot ACF and PACF for differenced avg_temp with short-term and long-term lags
plot_acf_pacf_diff(diff_avg_temp, "Differenced avg_temp", short_lag = 14, long_lag = 730)

# Plot ACF and PACF for differenced global_irradiation with short-term and long-term lags
plot_acf_pacf_diff(diff_global_irradiation, "Differenced global_irradiation", short_lag = 14, long_lag = 730)

# Plot ACF and PACF for differenced VOLUM_KWH with short-term and long-term lags
plot_acf_pacf_diff(diff_VOLUM_KWH, "Differenced VOLUM_KWH", short_lag = 14, long_lag = 730)
  

```

### Interpreting the differenced series

Compare these plots with the original ACF/PACF. Successful annual differencing should reduce peaks near lag 365 and its multiples; remaining short-lag structure can guide non-seasonal AR and MA terms. Avoid additional differencing unless the plots and stationarity diagnostics indicate that it is necessary.

### STL decomposition and smoothing

STL separates each daily series into seasonal, trend, and remainder components. A 365-day seasonal window models annual variation, a 547-day trend window smooths changes over roughly 18 months, and robust fitting reduces the influence of outliers. The deseasonalized series are then smoothed with a centered seven-day moving average.





```r

# Convert the original data to time series format with yearly seasonality (frequency = 365)
volum_kwh_ts <- ts(merged_data_clean$VOLUM_KWH, frequency = 365)
avg_temp_ts <- ts(merged_data_clean$avg_temp, frequency = 365)
global_irradiation_ts <- ts(merged_data_clean$global_irradiation, frequency = 365)

# Perform STL decomposition on the original data with sensible parameters
volum_kwh_stl <- stl(volum_kwh_ts, s.window = 365, t.window = 547, robust = TRUE)
avg_temp_stl <- stl(avg_temp_ts, s.window = 365, t.window = 547, robust = TRUE)
global_irradiation_stl <- stl(global_irradiation_ts, s.window = "periodic", t.window = 547, robust = TRUE)

# Extract seasonal components
volum_kwh_seasonal <- volum_kwh_stl$time.series[, "seasonal"]
avg_temp_seasonal <- avg_temp_stl$time.series[, "seasonal"]
global_irradiation_seasonal <- global_irradiation_stl$time.series[, "seasonal"]

# Create datasets without the seasonal component
volum_kwh_without_season <- volum_kwh_ts - volum_kwh_seasonal
avg_temp_without_season <- avg_temp_ts - avg_temp_seasonal
global_irradiation_without_season <- global_irradiation_ts - global_irradiation_seasonal

# Apply 7-day centered moving average to the deseasonalized data
volum_kwh_smoothed <- rollapply(volum_kwh_without_season, width = 7, FUN = mean, align = "center", fill = NA)
avg_temp_smoothed <- rollapply(avg_temp_without_season, width = 7, FUN = mean, align = "center", fill = NA)
global_irradiation_smoothed <- rollapply(global_irradiation_without_season, width = 7, FUN = mean, align = "center", fill = NA)

# Plot STL decomposition for each variable
plot(volum_kwh_stl, main = "STL Decomposition of VOLUM_KWH")
plot(avg_temp_stl, main = "STL Decomposition of avg_temp")
plot(global_irradiation_stl, main = "STL Decomposition of global_irradiation")


# Keep derived data in memory. Call write.csv() explicitly if an export is needed.
seasonal_components_df <- data.frame(
  Date = merged_data_clean$Date,
  volum_kwh_seasonal,
  avg_temp_seasonal,
  global_irradiation_seasonal
)
deseasonalized_smoothed_df <- data.frame(
  volum_kwh_smoothed,
  avg_temp_smoothed,
  global_irradiation_smoothed
)

```
```r

# Fill edge values introduced by the centered moving average with the nearest
# available smoothed value.
fill_edges <- function(x) {
  na.locf(na.locf(x, na.rm = FALSE), fromLast = TRUE, na.rm = FALSE)
}

deseasonalized_smoothed_df <- deseasonalized_smoothed_df %>%
  mutate(across(everything(), fill_edges))


# Check for any remaining NA values to ensure all are filled
print(sapply(deseasonalized_smoothed_df, function(x) sum(is.na(x))))
print(sapply(seasonal_components_df, function(x) sum(is.na(x))))

```



## Test lagged predictive relationships

A Granger test compares a restricted autoregression of a response against an unrestricted model that also contains lagged values of a candidate predictor. The null hypothesis is that those additional lags provide no predictive information. Rejection indicates incremental predictive value under the selected lag order, not physical causation. The same comparisons are run on the original and deseasonalized series.


```r

library(lmtest)

# Test up to 24 daily lags.
max_lag <- 24

# Apply the Granger causality test to the original data
granger_test_original_volum_kwh_avg_temp <- grangertest(
  merged_data_clean$VOLUM_KWH ~ merged_data_clean$avg_temp,
  order = max_lag
)
granger_test_original_volum_kwh_global_irradiation <- grangertest(
  merged_data_clean$VOLUM_KWH ~ merged_data_clean$global_irradiation,
  order = max_lag
)
granger_test_original_avg_temp_global_irradiation <- grangertest(
  merged_data_clean$avg_temp ~ merged_data_clean$global_irradiation,
  order = max_lag
)


# Print results for original data
print("Granger Causality Test on Original Data:")
print(granger_test_original_volum_kwh_avg_temp)
print(granger_test_original_volum_kwh_global_irradiation)
print(granger_test_original_avg_temp_global_irradiation)

# Apply the Granger causality test to the deseasonalized data
granger_test_deseasonalized_volum_kwh_avg_temp <- grangertest(
  deseasonalized_smoothed_df$volum_kwh_smoothed ~
    deseasonalized_smoothed_df$avg_temp_smoothed,
  order = max_lag
)
granger_test_deseasonalized_volum_kwh_global_irradiation <- grangertest(
  deseasonalized_smoothed_df$volum_kwh_smoothed ~
    deseasonalized_smoothed_df$global_irradiation_smoothed,
  order = max_lag
)
granger_test_deseasonalized_avg_temp_global_irradiation <- grangertest(
  deseasonalized_smoothed_df$avg_temp_smoothed ~
    deseasonalized_smoothed_df$global_irradiation_smoothed,
  order = max_lag
)

# Print results for deseasonalized data
print("Granger Causality Test on Deseasonalized Data:")
print(granger_test_deseasonalized_volum_kwh_avg_temp)
print(granger_test_deseasonalized_volum_kwh_global_irradiation)
print(granger_test_deseasonalized_avg_temp_global_irradiation)



```

### Interpret the Granger tests

For each test, reject the null only when the reported p-value is below the chosen significance level. A rejection means that lagged values of the predictor add predictive information for the response under the fitted lag specification; it does not establish a physical causal effect. Recheck conclusions after changing the lag order, sample period, or preprocessing.

## Evaluate univariate ETS forecasts

The following rolling-origin evaluation always trains on observations that precede the validation window. This chronological design avoids the future-to-past leakage that random cross-validation would introduce. ETS models provide a compact baseline for private consumption, temperature, and irradiation; weekly seasonality is used because the consumption series contains too few complete annual cycles for stable yearly estimation.


```r

# Rolling-window cross-validation preserves time order. The data are daily, so
# h, init_fold, and recent_window are expressed in days.
time_series_cv <- function(data, model_fct, init_fold, h, recent_window = Inf) {
  if (length(data) < init_fold + h) {
    stop("The series is too short for the requested initial fold and horizon.")
  }

  fold_starts <- seq(init_fold + 1, length(data) - h + 1, by = h)
  residuals <- lapply(fold_starts, function(fold_start) {
    train_start <- max(1, fold_start - min(recent_window, fold_start - 1))
    train <- ts(
      data[train_start:(fold_start - 1)],
      frequency = frequency(data)
    )
    test <- data[fold_start:(fold_start + h - 1)]
    predicted <- forecast(model_fct(train), h = h)$mean
    as.numeric(test - predicted)
  })

  residuals <- unlist(residuals)
  list(residuals = residuals, rmse = sqrt(mean(residuals^2)))
}

# Weekly seasonality is estimable by ETS for these daily series. Annual structure
# is handled separately by the STL template above.
train_model <- function(train_data) ets(train_data, model = "ZZZ")
volum_kwh_ts <- ts(merged_data_clean$VOLUM_KWH, frequency = 7)
temp_ts <- ts(merged_data_clean$avg_temp, frequency = 7)
irradiation_ts <- ts(merged_data_clean$global_irradiation, frequency = 7)

cv_settings <- list(init_fold = 365, h = 14, recent_window = 3 * 365)
volum_kwh_cv <- do.call(time_series_cv, c(list(volum_kwh_ts, train_model), cv_settings))
temp_cv <- do.call(time_series_cv, c(list(temp_ts, train_model), cv_settings))
irradiation_cv <- do.call(time_series_cv, c(list(irradiation_ts, train_model), cv_settings))

cat("CV RMSE — private electricity consumption:", volum_kwh_cv$rmse, "\n")
cat("CV RMSE — air temperature:", temp_cv$rmse, "\n")
cat("CV RMSE — global irradiation:", irradiation_cv$rmse, "\n")

forecast_horizon <- 14
volum_kwh_forecast <- forecast(train_model(volum_kwh_ts), h = forecast_horizon)
temp_forecast <- forecast(train_model(temp_ts), h = forecast_horizon)
irradiation_forecast <- forecast(train_model(irradiation_ts), h = forecast_horizon)

plot(volum_kwh_forecast, main = "14-Day Forecast: Private Electricity Consumption")
plot(temp_forecast, main = "14-Day Forecast: Air Temperature")
plot(irradiation_forecast, main = "14-Day Forecast: Global Irradiation")




volum_kwh_residuals <- volum_kwh_cv$residuals
temp_residuals <- temp_cv$residuals
irradiation_residuals <- irradiation_cv$residuals

# shapiro.test() accepts at most 5,000 values.
shapiro_volum_kwh <- shapiro.test(head(volum_kwh_residuals, 5000))
shapiro_temp_all <- shapiro.test(head(temp_residuals, 5000))
shapiro_irradiation_all <- shapiro.test(head(irradiation_residuals, 5000))

# Print Shapiro-Wilk test results for normality of residuals
cat("Shapiro-Wilk Test for Consumption Residuals: W =", shapiro_volum_kwh$statistic, ", p-value =", shapiro_volum_kwh$p.value, "\n")
cat("Shapiro-Wilk Test for Temperature Residuals (2017-2024): W =", shapiro_temp_all$statistic, ", p-value =", shapiro_temp_all$p.value, "\n")
cat("Shapiro-Wilk Test for Irradiation Residuals (2017-2024): W =", shapiro_irradiation_all$statistic, ", p-value =", shapiro_irradiation_all$p.value, "\n")


```

### Validation considerations

The Shapiro–Wilk tests are diagnostics rather than model-selection rules; inspect residual ACF and variance stability as well. With only a few annual cycles, unusual years can strongly influence estimated yearly seasonality and widen long-horizon uncertainty. The rolling validation window keeps the comparison realistic and limits the influence of remote history.

## Evaluate ARIMAX with a holdout period

The example below evaluates the final 90 observed days as a holdout period. A genuine future ARIMAX forecast requires forecasts or scenarios for temperature and irradiation; reusing the first 90 historical covariate rows would leak irrelevant information.

```r

arimax_holdout <- function(response, xreg, h = 90) {
  stopifnot(length(response) == nrow(xreg), length(response) > h)
  train_index <- seq_len(length(response) - h)
  test_index <- (length(response) - h + 1):length(response)

  model <- auto.arima(
    response[train_index],
    xreg = xreg[train_index, , drop = FALSE],
    seasonal = FALSE
  )
  predicted <- forecast(
    model,
    xreg = xreg[test_index, , drop = FALSE],
    h = h
  )

  list(
    model = model,
    forecast = predicted,
    actual = response[test_index],
    rmse = sqrt(mean((as.numeric(predicted$mean) - response[test_index])^2))
  )
}

raw_response <- ts(merged_data_clean$VOLUM_KWH, frequency = 365)
raw_xreg <- as.matrix(
  merged_data_clean[, c("avg_temp", "global_irradiation")]
)
raw_arimax <- arimax_holdout(raw_response, raw_xreg)

deseasonalized_response <- ts(
  deseasonalized_smoothed_df$volum_kwh_smoothed,
  frequency = 365
)
deseasonalized_xreg <- as.matrix(
  deseasonalized_smoothed_df[, c("avg_temp_smoothed", "global_irradiation_smoothed")]
)
deseasonalized_arimax <- arimax_holdout(
  deseasonalized_response,
  deseasonalized_xreg
)

cat("Raw-series ARIMAX holdout RMSE:", raw_arimax$rmse, "\n")
cat("Deseasonalized ARIMAX holdout RMSE:", deseasonalized_arimax$rmse, "\n")
plot(raw_arimax$forecast, main = "ARIMAX: 90-Day Holdout Forecast")
plot(deseasonalized_arimax$forecast, main = "Deseasonalized ARIMAX: 90-Day Holdout Forecast")

```

### Exogenous variables and deseasonalization

Temperature and irradiation can add useful information beyond the response history, but they also increase model complexity. Operational forecasts require future weather forecasts or explicit scenarios; observed future covariates must never leak into model training.

STL preprocessing lets ARIMAX focus on variation left after the estimated annual component and weekly smoothing. Compare the two holdout RMSE values before choosing a representation: deseasonalization may simplify the dynamics, but it can also remove predictive seasonal information.
