# %% [markdown]
# # Time-Series Similarity, Clustering, and Segmentation
#
# This template presents two complementary classification workflows. The first compares German electricity-generation series with several dissimilarity measures, hierarchical clustering, and hidden Markov segmentation. The second represents daily consumption as 24-hour profiles and evaluates whether unsupervised clusters recover weekday and weekend patterns. Generated plots and fitted output are intentionally omitted.
#
# %%
# Load the necessary packages
library(readr)
library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(TSclust)
library(pheatmap)
library(zoo)
library(tidyr)               
library(dtw)        
library(caret)
library(depmixS4)

resolve_data_file <- function(filename) {
  candidates <- c(filename, file.path("CA4-Classification", filename))
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop("Could not find ", filename, ". Run from the repository root or CA4-Classification/.")
  }

  existing[[1]]
}


# %% [markdown]
#
#
# ## Compare electricity-generation series
#
# ### Load and prepare generation data
#
# The source export contains quarter-hourly generation by production type. The loader accepts either a clean ENTSO-E export or the local course copy with a truncated duplicate prefix, then parses Berlin local time, removes unavailable source columns, and ranks the remaining sources by total production.
#
# %%

# The local course copy contains a truncated export followed by a clean
# full-year export. A fresh ENTSO-E download contains only the clean export.
file_path <- resolve_data_file(
  "Actual Generation per Production Type_202401010000-202501010000.csv"
)
generation_columns <- c(
  "Area", "MTU", "Biomass", "Fossil_Brown_Coal", "Fossil_Coal_Derived_Gas", "Fossil_Gas",
  "Fossil_Hard_Coal", "Fossil_Oil", "Fossil_Oil_Shale", "Fossil_Peat",
  "Geothermal", "Hydro_Pumped_Storage", "Hydro_Pumped_Storage_Consumption",
  "Hydro_Run_of_River_and_poundage", "Hydro_Water_Reservoir", "Marine", "Nuclear", "Other",
  "Other_Renewable", "Solar", "Waste", "Wind_Offshore", "Wind_Onshore"
)

read_generation_export <- function(file_path, column_names) {
  raw_lines <- readLines(file_path, warn = FALSE)
  header_pattern <- '"Area","MTU","Biomass'
  header_lines <- grep(header_pattern, raw_lines, fixed = TRUE)

  if (length(header_lines) == 0) {
    stop("Could not locate the ENTSO-E generation header in ", file_path)
  }

  last_header_line <- tail(header_lines, 1)
  has_clean_header <- startsWith(raw_lines[[last_header_line]], header_pattern)
  skip_lines <- if (has_clean_header) last_header_line - 1 else last_header_line

  generation_data <- read.csv(
    file_path,
    skip = skip_lines,
    header = has_clean_header,
    na.strings = c("n/e", "-", ""),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(generation_data) != length(column_names)) {
    stop("Unexpected number of columns in the ENTSO-E generation export.")
  }

  names(generation_data) <- column_names
  generation_data
}

generation_data <- read_generation_export(file_path, generation_columns)

# Change MTU to appropriate format
generation_data <- generation_data %>%
  mutate(
    MTU = str_split_fixed(MTU, " - ", 2)[, 1],  # Extract the start time
    MTU = dmy_hm(MTU, tz = "Europe/Berlin", quiet = TRUE)
  ) %>%
  # The source includes four nonexistent local times during the spring DST jump.
  filter(!is.na(MTU))

# Drop production types that are unavailable throughout this export.
generation_data <- generation_data %>%
  select(-c("Fossil_Coal_Derived_Gas", "Fossil_Oil_Shale", "Nuclear", "Fossil_Peat", "Marine"))




# Convert numerical columns to their appropriate format
generation_data <- generation_data %>%
  mutate(across(-c(Area, MTU), ~ as.numeric(.)))


# Exclude non-production-related columns
production_totals <- generation_data %>%
  summarise(across(-c(Area, MTU), ~ sum(.x, na.rm = TRUE)))


# Transform the summary into a long format for ranking
production_totals_long <- production_totals %>%
  pivot_longer(
    cols = everything(),
    names_to = "Source",
    values_to = "Total_Production"
  )


# Rank sources by total production and keep the top six.
top_sources <- production_totals_long %>%
  arrange(desc(Total_Production)) %>%
  slice_head(n = 6)


# Filter the dataset to keep only the top six sources
top_sources_list <- top_sources$Source
top_generation_data <- generation_data %>%
  select(all_of(c("Area", "MTU", top_sources_list)))

# Display the top contributors
print(top_sources)




# %% [markdown]
#
#
# ### Explore leading generation sources
#
# January provides a compact window for comparing the six largest sources. The faceted plot shows their original scales; the distance matrices then contrast magnitude, co-movement, autoregressive structure, and locally time-warped shape.
#
#
# %%


# Filter the data for January 2024
january_data <- top_generation_data %>%
  filter(
    MTU >= ymd_hms("2024-01-01 00:00:00", tz = "Europe/Berlin"),
    MTU < ymd_hms("2024-02-01 00:00:00", tz = "Europe/Berlin")
  ) %>%
  arrange(MTU) %>%
  mutate(across(all_of(top_sources_list), ~ na.approx(.x, na.rm = FALSE, rule = 2)))

# Reshape data for visualization
january_data_long <- january_data %>%
  pivot_longer(cols = -c(Area, MTU), names_to = "Source", values_to = "Production")

# Create time series plot with a grid layout
ggplot(january_data_long, aes(x = MTU, y = Production, color = Source)) +
  geom_line() +
  facet_wrap(~ Source, scales = "free_y", ncol = 2) + 
  theme_minimal() +
  labs(
    title = "Electricity Production by Top 6 Sources (January 2024)",
    x = "Date and Time",
    y = "Production (MW)"
  ) +
  theme(legend.position = "none")








# Select relevant columns 
generation_series <- january_data %>% select(-c(Area, MTU))

# Transpose the matrix so rows become production types
transposed_columns <- t(as.matrix(generation_series))


# Euclidean Distance
euclidean_dist <- diss(transposed_columns, "EUCL")

# Pearson Correlation Dissimilarity (1 - Correlation)
pearson_dissimilarity <- diss(transposed_columns, "COR")

# AR.PIC 
ar_pic_dist <- diss(transposed_columns, "AR.PIC")

# Dynamic Time Warping (DTW) Dissimilarity
dtw_dist <- diss(transposed_columns, "DTW")



# Euclidean Distance Heatmap
pheatmap(as.matrix(euclidean_dist),
         main = "Euclidean Distance Between Production Types",
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# Pearson Correlation Dissimilarity Heatmap
pheatmap(as.matrix(pearson_dissimilarity),
         main = "Pearson Correlation Dissimilarity Between Production Types",
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# AR.PIC Dissimilarity Heatmap
pheatmap(as.matrix(ar_pic_dist),
         main = "AR.PIC Dissimilarity Between Production Types",
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# DTW Dissimilarity Heatmap
pheatmap(as.matrix(dtw_dist),
         main = "DTW Dissimilarity Between Production Types",
         cluster_rows = TRUE,
         cluster_cols = TRUE)





# Compute Autocorrelations for Each Series
autocorrelations <- lapply(generation_series, function(series) {
  acf(series, plot = FALSE)  
})

# Convert Autocorrelations into a Data Frame for Visualization
autocorrelation_data <- do.call(rbind, lapply(seq_along(autocorrelations), function(i) {
  data.frame(
    Lag = autocorrelations[[i]]$lag,  
    ACF = autocorrelations[[i]]$acf,      
    Source = names(generation_series)[i]
  )
}))


# Number of observations used for approximate ACF confidence bounds.
n_observations <- nrow(generation_series)


# Visualize Autocorrelations in a Grid of Plots
ggplot(autocorrelation_data, aes(x = Lag, y = ACF)) +
  geom_col(position = "dodge", fill = "steelblue") +
  geom_hline(
    yintercept = c(-1.96 / sqrt(n_observations), 1.96 / sqrt(n_observations)),
    linetype = "dashed",
    color = "red"
  ) +
  facet_wrap(~ Source, ncol = 2) +  # Create a grid of plots
  theme_minimal() +
  labs(
    title = "Autocorrelation For Top 6 Electricity Sources",
    x = "Lag",
    y = "Autocorrelation"
  )



# %% [markdown]
# ### Summary of the exploratory analysis
#
# Use Euclidean distance to compare absolute production magnitudes, correlation dissimilarity to compare scale-free co-movement, DTW to compare shapes with local timing shifts, and AR.PIC to compare autoregressive dynamics. The ACF panels show persistence within each source. Run the paired script to identify the source-specific groupings; saved results are intentionally omitted.
#
# Euclidean distance and Pearson correlation dissimilarity are carried forward to show how a scale-sensitive and scale-free definition of similarity produce different hierarchies.
#
#
#
# ### Hierarchical clustering
#
# %%

# Perform hierarchical clustering 
hc_euclidean <- hclust(dist(t(as.matrix(generation_series))), method = "ward.D2")
hc_pearson <- hclust(as.dist(1 - cor(generation_series)), method = "average")


# Plot dendrograms
plot(hc_euclidean, 
     main = "Hierarchical Clustering (Euclidean Distance)", 
     xlab = "Electricity Sources", 
     sub = "", 
     hang = -1)

plot(hc_pearson, 
     main = "Hierarchical Clustering (Pearson Correlation Dissimilarity)", 
     xlab = "Electricity Sources", 
     sub = "", 
     hang = -1)



# Assign clusters to columns
clusters_euclidean <- cutree(hc_euclidean, k = 3)  
clusters_pearson <- cutree(hc_pearson, k = 3) 


# Map clusters to column names 
clusters_euclidean_df <- data.frame(
  Source = colnames(generation_series),
  Cluster = as.factor(clusters_euclidean)
)
clusters_pearson_df <- data.frame(
  Source = colnames(generation_series),
  Cluster = as.factor(clusters_pearson)
)


# Reshape the January data for visualization
january_data_long <- january_data %>%
  pivot_longer(cols = -c(Area, MTU), names_to = "Source", values_to = "Production")


# Map clusters back to the reshaped data
january_data_long_euclidean <- january_data_long %>%
  left_join(clusters_euclidean_df, by = "Source") %>%
  rename(Euclidean_Cluster = Cluster)

january_data_long_pearson <- january_data_long %>%
  left_join(clusters_pearson_df, by = "Source") %>%
  rename(Pearson_Cluster = Cluster)



# Plotting time series
ggplot(january_data_long_euclidean, aes(x = MTU, y = Production, color = Euclidean_Cluster, group = Source)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Source, ncol = 2, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Time Series Colored by Euclidean Clusters",
    x = "Date and Time",
    y = "Production (MW)"
  ) +
  theme(legend.position = "bottom")

ggplot(january_data_long_pearson, aes(x = MTU, y = Production, color = Pearson_Cluster, group = Source)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Source, ncol = 2, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Time Series Colored by Pearson Clusters",
    x = "Date and Time",
    y = "Production (MW)"
  ) +
  theme(legend.position = "bottom")





# %% [markdown]
# ### Interpretation of the groupings
#
# Compare the dendrograms to see how scale-sensitive Euclidean distance and scale-free correlation dissimilarity change the groupings. Interpret cluster membership only after regenerating the plots, because the document stores no fitted output.
#
# ### Hidden Markov segmentation
#
# A separate two-state Gaussian hidden Markov model divides each generation series into latent operating regimes. The states are data-driven labels rather than predefined low/high categories, so interpret them from the regenerated state-colored plots.
#
# %%

# Create a Function to fit a two-state HMM 
fit_hmm <- function(time_series) {
  set.seed(123)
  model <- depmix(response = time_series ~ 1,
                  family = gaussian(),
                  nstates = 2,
                  data = data.frame(time_series = time_series))
  
  fit <- fit(model)
  return(fit)
}

# Apply the HMM fitting function 
hmm_models <- lapply(generation_series, fit_hmm)

# Function to extract the most likely states from an HMM
extract_states <- function(model) {
  posterior(model, type = "viterbi")$state
}

# Extract states for each model
state_sequences <- lapply(hmm_models, extract_states)

# Combine state sequences into a data frame
state_data <- data.frame(MTU = january_data$MTU, do.call(cbind, state_sequences))
colnames(state_data)[-1] <- colnames(generation_series)





# %%

# Reshape state data for visualization
state_data_long <- state_data %>%
  pivot_longer(cols = -MTU, names_to = "Source", values_to = "State") %>%
  left_join(january_data_long, by = c("MTU", "Source"))


# Plot time series with state segmentation
ggplot(state_data_long, aes(x = MTU, y = Production, color = factor(State), group = Source)) +
  geom_line() +
  facet_wrap(~ Source, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("blue", "red"), labels = c("State 1", "State 2")) +
  theme_minimal() +
  labs(
    title = "Time Series Colored by Hidden Markov Model States",
    x = "Date and Time",
    y = "Production (MW)",
    color = "State"
  )


# %% [markdown]
#
# ### How the generation methods differ
#
# The methods answer complementary questions: distance measures compare whole series, hierarchical clustering summarizes between-series similarity, autocorrelation describes within-series persistence, and hidden Markov models segment each series into latent regimes. Euclidean distance is sensitive to scale, correlation focuses on co-movement, and DTW tolerates local timing shifts. Use the regenerated plots and fitted state sequences to make source-specific claims.
#
# ## Classify weekday and weekend consumption
#
# ### Load and summarize consumption
#
# The consumption data contain hourly readings for private, business, and industrial groups. After converting timestamps to Oslo local time, the code selects the longest contiguous date range and summarizes mean daily consumption by weekday. Averaging daily totals avoids confounding weekday patterns with unequal numbers of Mondays, Tuesdays, and other weekdays in the sample.
#
# %%

# Import the dataset. read_csv2() handles decimal commas and ISO-8601 offsets.
file_path <- resolve_data_file("consumption_per_group_aas_hour.csv")
consumption_data <- read_csv2(
  file_path,
  locale = locale(tz = "Europe/Oslo"),
  show_col_types = FALSE
)

# The course copy is already limited to Ås; retain the same scope if a full
# Elhub municipality download is supplied instead.
if ("KOMMUNE" %in% names(consumption_data)) {
  consumption_data <- consumption_data %>%
    filter(KOMMUNE == "Ås")
}

consumption_data <- consumption_data %>%
  mutate(
    STARTTID = with_tz(STARTTID, "Europe/Oslo"),
    SLUTTID = with_tz(SLUTTID, "Europe/Oslo")
  )

# Missing Values
missing_values_consum <- sapply(consumption_data, function(x) sum(is.na(x)))
print("Missing Values Per Column:")
print(missing_values_consum)


# Retain the longest contiguous suffix instead of relying on a hard-coded date.
consumption_dates <- sort(unique(as.Date(consumption_data$STARTTID)))
gap_indices <- which(diff(consumption_dates) > 1)
contiguous_start <- if (length(gap_indices) > 0) {
  consumption_dates[tail(gap_indices, 1) + 1]
} else {
  consumption_dates[1]
}

consumption_data <- consumption_data %>%
  filter(as.Date(STARTTID) >= contiguous_start)





# Create local calendar fields and daily totals.
daily_consumption <- consumption_data %>%
  mutate(
    Date = as.Date(STARTTID),
    Weekday = wday(
      STARTTID,
      label = TRUE,
      abbr = FALSE,
      week_start = 1,
      locale = "C"
    )
  ) %>%
  group_by(Date, FORBRUKSGRUPPE, Weekday) %>%
  summarise(Daily_Consumption = sum(VOLUM_KWH, na.rm = TRUE), .groups = "drop")

# Average comparable daily totals for each weekday and consumer group.
weekday_summary <- daily_consumption %>%
  group_by(FORBRUKSGRUPPE, Weekday) %>%
  summarize(
    Mean_Daily_Consumption = mean(Daily_Consumption),
    .groups = "drop"
  )

# Plot for Industry
industry_data <- weekday_summary %>%
  filter(FORBRUKSGRUPPE == "Industri")

industry_plot <- ggplot(
  industry_data,
  aes(x = Weekday, y = Mean_Daily_Consumption, group = 1)
) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Weekly Energy Consumption Patterns - Industry",
    x = "Weekday",
    y = "Mean Daily Consumption (kWh)"
  ) +
  theme_minimal()

# Plot for Business
business_data <- weekday_summary %>%
  filter(FORBRUKSGRUPPE == "Forretning")

business_plot <- ggplot(
  business_data,
  aes(x = Weekday, y = Mean_Daily_Consumption, group = 1)
) +
  geom_line(color = "green", linewidth = 1.2) +
  geom_point(color = "purple", size = 2) +
  labs(
    title = "Weekly Energy Consumption Patterns - Business",
    x = "Weekday",
    y = "Mean Daily Consumption (kWh)"
  ) +
  theme_minimal()

# Plot for Private
private_data <- weekday_summary %>%
  filter(FORBRUKSGRUPPE == "Privat")

private_plot <- ggplot(
  private_data,
  aes(x = Weekday, y = Mean_Daily_Consumption, group = 1)
) +
  geom_line(color = "orange", linewidth = 1.2) +
  geom_point(color = "brown", size = 2) +
  labs(
    title = "Weekly Energy Consumption Patterns - Private",
    x = "Weekday",
    y = "Mean Daily Consumption (kWh)"
  ) +
  theme_minimal()

# Display the plots
print(industry_plot)
print(business_plot)
print(private_plot)




# %% [markdown]
#
# ### Compare daily load shapes
#
# Each consumer group is reshaped into a matrix whose columns are complete days and whose rows are local clock hours. Euclidean distance compares profiles hour by hour, while dynamic time warping (DTW) tolerates small shifts in the timing of peaks. The eight-week window avoids daylight-saving days, which would otherwise require explicit handling of 23- and 25-hour profiles.
#
# %%



# Limit data to a specific 8-week period
start_date <- as.Date("2024-01-01") 
end_date <- start_date + 56          

data_8weeks <- consumption_data %>%
  filter(as.Date(STARTTID) >= start_date & as.Date(STARTTID) < end_date)


# Function to transform data into a 24 x n matrix for each consumer group
transform_to_matrix <- function(input_data, group_name) {
  group_data <- input_data %>%
    filter(FORBRUKSGRUPPE == group_name) %>%
    mutate(
      Hour = format(STARTTID, "%H"), 
      Date = as.Date(STARTTID)
    ) %>%
    group_by(Date, Hour) %>%
    summarize(
      Total_Consumption = sum(VOLUM_KWH, na.rm = TRUE), 
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Date,
      values_from = Total_Consumption,
      values_fill = NA_real_
    ) %>%
    arrange(as.numeric(Hour))
  
  # Convert to matrix and set row names as hours
  matrix_data <- as.matrix(group_data[ , -1])
  rownames(matrix_data) <- group_data$Hour
  complete_days <- colSums(is.na(matrix_data)) == 0
  matrix_data[, complete_days, drop = FALSE]
}



# Function to compute pairwise distance matrix 
compute_daily_distance_matrix <- function(data_matrix, method) {
  method <- match.arg(method, c("euclidean", "dtw"))

  if (method == "euclidean") {
    return(as.matrix(dist(t(data_matrix), method = "euclidean")))
  }
  
  n <- ncol(data_matrix)  # Number of days 
  distance_matrix <- matrix(0, nrow = n, ncol = n)
  colnames(distance_matrix) <- colnames(data_matrix)
  rownames(distance_matrix) <- colnames(data_matrix)
  
  # Loop through all pairs of days
  for (i in 1:n) {
    for (j in i:n) {  # Compute distances only for unique pairs
      dist_val <- dtw(
        data_matrix[, i],
        data_matrix[, j],
        distance.only = TRUE
      )$distance
      distance_matrix[i, j] <- dist_val
      distance_matrix[j, i] <- dist_val  
    }
  }
  return(distance_matrix)
}


smooth_matrix <- function(data_matrix, k = 3) {
  if (k < 1 || k > nrow(data_matrix)) {
    stop("k must be between 1 and the number of hourly rows.")
  }

  apply(data_matrix, 2, function(column) {
    zoo::rollmean(column, k = k, fill = "extend", align = "center")
  })
}

normalize <- function(x) {
  value_range <- range(x, na.rm = TRUE)
  if (diff(value_range) == 0) {
    return(rep(0, length(x)))
  }
  (x - value_range[1]) / diff(value_range)
}

# Function to visualize the distance matrix as a heatmap
visualize_distance_matrix <- function(distance_matrix, title, palette_function) {
  pheatmap(
    distance_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = title,
    color = palette_function(100),
    border_color = NA
  )
}


# Create 24 x n Matrices for Each Consumer Group 
private_matrix <- transform_to_matrix(data_8weeks, "Privat")
business_matrix <- transform_to_matrix(data_8weeks, "Forretning")
industry_matrix <- transform_to_matrix(data_8weeks, "Industri")


# Compute Distance Matrices 
private_euclidean <- compute_daily_distance_matrix(private_matrix, method = "euclidean")
private_dtw <- compute_daily_distance_matrix(private_matrix, method = "dtw")

business_euclidean <- compute_daily_distance_matrix(business_matrix, method = "euclidean")
business_dtw <- compute_daily_distance_matrix(business_matrix, method = "dtw")

industry_euclidean <- compute_daily_distance_matrix(industry_matrix, method = "euclidean")
industry_dtw <- compute_daily_distance_matrix(industry_matrix, method = "dtw")


# Visualize Distance Matrices
visualize_distance_matrix(private_euclidean, title = "Euclidean Distance (Private - 8 Weeks)",
                          palette_function = viridis::viridis)
visualize_distance_matrix(private_dtw, title = "DTW Distance (Private - 8 Weeks)",
                          palette_function = viridis::viridis)

visualize_distance_matrix(business_euclidean, title = "Euclidean Distance (Business - 8 Weeks)",
                          palette_function = viridis::magma)
visualize_distance_matrix(business_dtw, title = "DTW Distance (Business - 8 Weeks)",
                          palette_function = viridis::magma)

visualize_distance_matrix(industry_euclidean, title = "Euclidean Distance (Industry - 8 Weeks)",
                          palette_function = viridis::cividis)
visualize_distance_matrix(industry_dtw, title = "DTW Distance (Industry - 8 Weeks)",
                          palette_function = viridis::cividis)





# Function to generate and plot a dendrogram 
visualize_dendrogram <- function(distance_matrix, title) {
  hc <- hclust(as.dist(distance_matrix), method = "complete")

  plot(
    hc,
    main = title,
    xlab = "Days",
    sub = "",
    cex.main = 1.5,
    cex.lab = 1.2,
    cex.axis = 0.8
  )
}


# Dendogram Visualization 
visualize_dendrogram(private_euclidean, title = "Dendrogram (Private - Euclidean Distance)")
visualize_dendrogram(private_dtw, title = "Dendrogram (Private - DTW Distance)")

visualize_dendrogram(business_euclidean, title = "Dendrogram (Business - Euclidean Distance)")
visualize_dendrogram(business_dtw, title = "Dendrogram (Business - DTW Distance)")

visualize_dendrogram(industry_euclidean, title = "Dendrogram (Industry - Euclidean Distance)")
visualize_dendrogram(industry_dtw, title = "Dendrogram (Industry - DTW Distance)")






# Function to visualize the distance matrix as a heatmap (chronological)
visualize_chronological_heatmap <- function(distance_matrix, title, palette_function) {
  pheatmap(
    distance_matrix,
    cluster_rows = FALSE,  
    cluster_cols = FALSE,  
    main = title,
    color = palette_function(100),  
    border_color = NA
  )
}


# Private consumer group
visualize_chronological_heatmap(private_euclidean, 
                                title = "Euclidean Distance (Private - Chronological)", 
                                palette_function = viridis::viridis)
visualize_chronological_heatmap(private_dtw, 
                                title = "DTW Distance (Private - Chronological)", 
                                palette_function = viridis::viridis)

# Business consumer group
visualize_chronological_heatmap(business_euclidean, 
                                title = "Euclidean Distance (Business - Chronological)", 
                                palette_function = viridis::magma)
visualize_chronological_heatmap(business_dtw, 
                                title = "DTW Distance (Business - Chronological)", 
                                palette_function = viridis::magma)

# Industry consumer group
visualize_chronological_heatmap(industry_euclidean, 
                                title = "Euclidean Distance (Industry - Chronological)", 
                                palette_function = viridis::cividis)
visualize_chronological_heatmap(industry_dtw, 
                                title = "DTW Distance (Industry - Chronological)", 
                                palette_function = viridis::cividis)



# %% [markdown]
#
#
# ### Cluster and evaluate daily profiles
#
# The profiles are smoothed with a three-hour centered mean and normalized within each day so clustering emphasizes shape rather than absolute consumption. Cluster labels are arbitrary, so the evaluation checks both possible label mappings before calculating weekday/weekend accuracy. Ward's method is used for Euclidean distances; average linkage is used for DTW because Ward's variance criterion assumes Euclidean geometry.
#
# %%

# Smooth the matrices with a rolling mean (k = 3)
private_matrix <- smooth_matrix(private_matrix, k = 3)
business_matrix <- smooth_matrix(business_matrix, k = 3)
industry_matrix <- smooth_matrix(industry_matrix, k = 3)

# Normalize the matrices (scaling between 0 and 1)
private_matrix <- apply(private_matrix, 2, normalize)
business_matrix <- apply(business_matrix, 2, normalize)
industry_matrix <- apply(industry_matrix, 2, normalize)

# Compute distance matrices
private_euclidean <- compute_daily_distance_matrix(private_matrix, method = "euclidean")
business_euclidean <- compute_daily_distance_matrix(business_matrix, method = "euclidean")
industry_euclidean <- compute_daily_distance_matrix(industry_matrix, method = "euclidean")

private_dtw <- compute_daily_distance_matrix(private_matrix, method = "dtw")
business_dtw <- compute_daily_distance_matrix(business_matrix, method = "dtw")
industry_dtw <- compute_daily_distance_matrix(industry_matrix, method = "dtw")

# Perform hierarchical clustering
hclust_euclidean_priv <- hclust(as.dist(private_euclidean), method = "ward.D2")
hclust_dtw_priv <- hclust(as.dist(private_dtw), method = "average")

hclust_euclidean_bus <- hclust(as.dist(business_euclidean), method = "ward.D2")
hclust_dtw_bus <- hclust(as.dist(business_dtw), method = "average")

hclust_euclidean_ind <- hclust(as.dist(industry_euclidean), method = "ward.D2")
hclust_dtw_ind <- hclust(as.dist(industry_dtw), method = "average")

# Cut the dendrogram into two clusters
clusters_euclidean_priv <- cutree(hclust_euclidean_priv, k = 2)
clusters_dtw_priv <- cutree(hclust_dtw_priv, k = 2)

clusters_euclidean_bus <- cutree(hclust_euclidean_bus, k = 2)
clusters_dtw_bus <- cutree(hclust_dtw_bus, k = 2)

clusters_euclidean_ind <- cutree(hclust_euclidean_ind, k = 2)
clusters_dtw_ind <- cutree(hclust_dtw_ind, k = 2)

align_cluster_labels <- function(clusters, truth) {
  flipped <- ifelse(clusters == 1, 2, 1)
  if (mean(clusters == truth) >= mean(flipped == truth)) clusters else flipped
}

# Cluster IDs are arbitrary. Align the two labels to the weekday/weekend truth
# before computing accuracy; otherwise a perfectly inverted solution scores 0%.
create_confusion <- function(clusters, data_matrix) {
  truth <- ifelse(
    wday(as.Date(colnames(data_matrix))) %in% c(1, 7),
    2,
    1
  )
  aligned <- align_cluster_labels(clusters, truth)

  confusionMatrix(
    factor(aligned, levels = c(1, 2)),
    factor(truth, levels = c(1, 2))
  )
}


# Private sector
confusion_euclidean_priv <- create_confusion(clusters_euclidean_priv, private_matrix)
confusion_dtw_priv <- create_confusion(clusters_dtw_priv, private_matrix)

# Business sector
confusion_euclidean_bus <- create_confusion(clusters_euclidean_bus, business_matrix)
confusion_dtw_bus <- create_confusion(clusters_dtw_bus, business_matrix)

# Industry sector
confusion_euclidean_ind <- create_confusion(clusters_euclidean_ind, industry_matrix)
confusion_dtw_ind <- create_confusion(clusters_dtw_ind, industry_matrix)




# Print confusion matrices
print("Private - Euclidean")
print(confusion_euclidean_priv)

print("Private - DTW")
print(confusion_dtw_priv)

print("Business - Euclidean")
print(confusion_euclidean_bus)

print("Business - DTW")
print(confusion_dtw_bus)

print("Industry - Euclidean")
print(confusion_euclidean_ind)

print("Industry - DTW")
print(confusion_dtw_ind)





# Visualize Clustering Results
metrics <- data.frame(
  Group = rep(c("Private", "Business", "Industry"), each = 2),
  Method = rep(c("Euclidean", "DTW"), 3),
  Accuracy = c(
    confusion_euclidean_priv$overall["Accuracy"], confusion_dtw_priv$overall["Accuracy"],
    confusion_euclidean_bus$overall["Accuracy"], confusion_dtw_bus$overall["Accuracy"],
    confusion_euclidean_ind$overall["Accuracy"], confusion_dtw_ind$overall["Accuracy"]
  )
)

ggplot(metrics, aes(x = Group, y = Accuracy, fill = Method)) +
  geom_col(position = "dodge") +
  labs(title = "Clustering Accuracy by Consumer Group and Distance Measure",
       x = "Consumer Group", y = "Accuracy") +
  theme_minimal()



# %% [markdown]
#
# #### Interpret clustering metrics
# The code prints confusion matrices and plots accuracy after aligning arbitrary cluster IDs to the weekday/weekend labels. Compare Euclidean and DTW within each consumer group rather than treating the raw cluster numbers as class labels. Euclidean distance emphasizes point-by-point magnitude differences, whereas DTW permits local temporal alignment and may be more useful when peaks occur at slightly different hours.
#
# Because generated outputs are intentionally omitted, rerun the analysis before drawing conclusions about which method or consumer group performs best.
#
#
#
# ### Diagnose misclassified days
#
# The diagnostic below focuses on the industry DTW clustering. Comparing misclassified and correctly classified daily curves can reveal where temporal alignment is insufficient. The eight-week window begins on 1 January, so public holidays must be considered explicitly.
#  
# %%

# Ground truth (weekday = 1, weekend = 2)
ground_truth <- ifelse(wday(as.Date(colnames(industry_matrix))) %in% c(1, 7), 2, 1)

# Predictions for the industry group (DTW)
predicted_clusters <- align_cluster_labels(clusters_dtw_ind, ground_truth)

# Find misclassified and correctly classified days
misclassified_days <- colnames(industry_matrix)[predicted_clusters != ground_truth]
correctly_classified_days <- colnames(industry_matrix)[predicted_clusters == ground_truth]

# Print misclassified days
print("Misclassified Days:")
print(misclassified_days)
print("Correctly Classified days:")
print(correctly_classified_days)



# Extract profiles for diagnostic comparison.
misclassified_data <- industry_matrix[
  , colnames(industry_matrix) %in% misclassified_days, drop = FALSE
]
correctly_classified_data <- industry_matrix[
  , colnames(industry_matrix) %in% correctly_classified_days, drop = FALSE
]

if (length(misclassified_days) > 0 && length(correctly_classified_days) > 0) {
  matplot(
    0:23, misclassified_data, type = "l", col = "red", lty = 1, lwd = 2,
    xlab = "Hour", ylab = "Normalized consumption",
    main = "Misclassified and Correctly Classified Industry Days",
    ylim = range(industry_matrix)
  )
  matlines(
    0:23, correctly_classified_data, col = "blue", lty = 2, lwd = 1
  )
  legend(
    "topleft",
    legend = c("Misclassified", "Correctly classified"),
    col = c("red", "blue"),
    lty = c(1, 2),
    lwd = c(2, 1)
  )
} else {
  message("Both correct and incorrect classifications are required for profile comparison.")
}




# Norwegian public holidays for 2024
norwegian_holidays <- as.Date(c(
  "2024-01-01",  # New Year's Day
  "2024-03-28",  # Maundy Thursday
  "2024-03-29",  # Good Friday
  "2024-03-31",  # Easter Sunday
  "2024-04-01",  # Easter Monday
  "2024-05-01",  # Labour Day
  "2024-05-17",  # Constitution Day
  "2024-05-09",  # Ascension Day
  "2024-05-19",  # Whit Sunday
  "2024-05-20",  # Whit Monday
  "2024-12-25",  # Christmas Day
  "2024-12-26"   # Boxing Day
))


# Convert misclassified days to Date format
misclassified_days_as_date <- as.Date(misclassified_days)

# Identify holidays among misclassified days
special_days <- misclassified_days_as_date[misclassified_days_as_date %in% norwegian_holidays]

print("Special Days (Holidays) among Misclassified Days:")
print(special_days)






# %%

if (length(misclassified_days) > 0 && length(correctly_classified_days) > 0) {
  profile_summary <- data.frame(
    Hour = 0:23,
    AvgMisclassified = rowMeans(misclassified_data, na.rm = TRUE),
    SDMisclassified = apply(misclassified_data, 1, sd, na.rm = TRUE),
    AvgCorrectlyClassified = rowMeans(correctly_classified_data, na.rm = TRUE),
    SDCorrectlyClassified = apply(correctly_classified_data, 1, sd, na.rm = TRUE)
  )

  ggplot(profile_summary, aes(x = Hour)) +
    geom_ribbon(
      aes(
        ymin = AvgMisclassified - SDMisclassified,
        ymax = AvgMisclassified + SDMisclassified
      ),
      fill = "red",
      alpha = 0.2
    ) +
    geom_ribbon(
      aes(
        ymin = AvgCorrectlyClassified - SDCorrectlyClassified,
        ymax = AvgCorrectlyClassified + SDCorrectlyClassified
      ),
      fill = "blue",
      alpha = 0.2
    ) +
    geom_line(aes(y = AvgMisclassified), color = "red", linewidth = 1) +
    geom_line(aes(y = AvgCorrectlyClassified), color = "blue", linewidth = 1) +
    labs(
      title = "Average Industry Profiles with Hourly Variability",
      x = "Hour",
      y = "Normalized consumption"
    ) +
    theme_minimal()
}


# %% [markdown]
#
# The profile comparison helps distinguish shape overlap from calendar effects. Cross-reference any errors with `special_days`: a public holiday can follow a weekend-like demand pattern even when its calendar label is a weekday. Do not infer a source-specific explanation until the diagnostics have been regenerated.
