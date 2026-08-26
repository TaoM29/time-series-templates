# Air-quality preprocessing and PCA

This template prepares hourly air-quality sensor measurements for exploration and principal component analysis (PCA). It compares complete-case PCA with a version fitted after linear interpolation, then reconstructs the measurements from a reduced set of components. Generated results and plots are intentionally omitted; run the paired `task1.R` script to reproduce them.

```r


# Air-quality preprocessing and PCA template
#
# Run this script from the repository root or from CA2-Preprocessing/.
library(tsibble)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(factoextra)
library(reshape2)
library(zoo)

resolve_data_file <- function(filename) {
  candidates <- c(filename, file.path("CA2-Preprocessing", filename))
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop("Could not find ", filename, ". Run from the repository root or CA2-Preprocessing/.")
  }

  existing[[1]]
}


```

## Load and explore hourly measurements

The UCI dataset uses `-200` for missing observations and stores dates and times separately. The following preparation combines them into a Rome-local timestamp, removes empty columns and incomplete calendar days, converts sensor fields to numeric values, and plots each measurement on its own scale.

```r

# Load the dataset
data <- read.table(
  resolve_data_file("AirQualityUCI.csv"),
  sep = ";",
  header = TRUE,
  na.strings = "-200"
)


# Remove columns that contain only NA values
data <- data %>%
  select(-X, -X.1)



# Combine Date and Time into a single DateTime column
data <- data %>%
  mutate(
    DateTime = parse_date_time(
      paste(Date, Time),
      orders = "dmy HMS",
      tz = "Europe/Rome"
    )
  ) %>%
  select(DateTime, everything(), -Date, -Time)  # Ensure Date and Time are removed



# Remove rows with NA in DateTime column
data <- data %>%
  filter(!is.na(DateTime))


# Convert to tsibble with DateTime as index
data_tsibble <- data %>%
  as_tsibble(index = DateTime)



# Identify complete days (days with exactly 24 hourly records)
data_tsibble <- data_tsibble %>%
  group_by(Date = as.Date(DateTime)) %>%  
  filter(n() == 24) %>%
  ungroup()


# Explicitly remove the Date column after grouping
data_tsibble <- data_tsibble %>%
  select(-Date)



# Convert variables to appropriate types
data_tsibble <- data_tsibble %>%
  mutate(
    CO.GT. = as.numeric(gsub(",", ".", as.character(CO.GT.))),
    PT08.S1.CO. = as.integer(PT08.S1.CO.),
    NMHC.GT. = as.integer(NMHC.GT.),
    C6H6.GT. = as.numeric(gsub(",", ".", as.character(C6H6.GT.))),
    PT08.S2.NMHC. = as.integer(PT08.S2.NMHC.),
    NOx.GT. = as.integer(ifelse(NOx.GT. == -200, NA, NOx.GT.)),  # Handle -200 as NA
    PT08.S3.NOx. = as.integer(PT08.S3.NOx.),
    NO2.GT. = as.integer(ifelse(NO2.GT. == -200, NA, NO2.GT.)),  # Handle -200 as NA
    PT08.S4.NO2. = as.integer(PT08.S4.NO2.),
    PT08.S5.O3. = as.integer(PT08.S5.O3.),
    T = as.numeric(gsub(",", ".", as.character(T))),
    RH = as.numeric(gsub(",", ".", as.character(RH))),
    AH = as.numeric(gsub(",", ".", as.character(AH)))
  )


# Gather the data to long format for plotting, excluding DateTime
data_long <- data_tsibble %>%
  pivot_longer(
    cols = -DateTime,  # Exclude DateTime
    names_to = "Variable", 
    values_to = "Value"
  )


# Plotting
ggplot(data_long, aes(x = DateTime, y = Value, color = Variable)) +
  geom_line() +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +  # Create multiple panels
  labs(title = "Air Quality Measurements Over Time",
       x = "DateTime",
       y = "Value") +
  theme_minimal() + 
  theme(legend.position = "none")



```

## Establish a complete-case PCA baseline

PCA requires a complete numeric matrix. This baseline keeps only rows without missing measurements, standardizes variables to comparable scales, and visualizes component variance, loadings, and scores.

```r

# Convert tsibble to a regular data frame
data_df <- as.data.frame(data_tsibble)


# Remove DateTime and prepare for PCA
data_for_pca <- data_df %>%
  select(-DateTime) 


# Remove any rows with NA values to ensure PCA can run
data_for_pca <- na.omit(data_for_pca)


# Perform PCA
pca_result <- prcomp(data_for_pca, center = TRUE, scale. = TRUE)


# Create Scree Plot with correctly ordered PC factors
screeplot_data <- data.frame(
  PC = factor(paste0("PC", 1:length(pca_result$sdev)), levels = paste0("PC", 1:length(pca_result$sdev))),
  Variance = pca_result$sdev^2
)

ggplot(screeplot_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(y = cumsum(Variance) / sum(Variance) * max(Variance)), color = "red") +
  labs(title = "Scree Plot", x = "Principal Components", y = "Variance") +
  theme_minimal()



# Biplot for PC1 and PC2
fviz_pca_biplot(pca_result, 
                repel = TRUE, 
                title = "PCA Biplot: PC1 vs PC2",
                col.var = "blue", 
                col.ind = adjustcolor("black", alpha.f = 0.3)) 


# Biplot for PC2 and PC3
fviz_pca_biplot(pca_result, axes = c(2, 3), 
                repel = TRUE, 
                title = "PCA Biplot: PC2 vs PC3",
                col.var = "blue", 
                col.ind = adjustcolor("black", alpha.f = 0.3)) 



# Plot the Scores for the Principal Components
scores <- as.data.frame(pca_result$x)
scores <- scores %>% mutate(Sample = row_number())  # Add a sample identifier


# Plot the scores for the first two principal components
ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "Scores Plot: PC1 vs PC2", x = "PC1", y = "PC2") +
  theme_minimal()


# Plot the scores for the second and third principal components
ggplot(scores, aes(x = PC2, y = PC3)) +
  geom_point() +
  labs(title = "Scores Plot: PC2 vs PC3", x = "PC2", y = "PC3") +
  theme_minimal()



```

## Diagnose and impute missing values

Correlated missingness can indicate sensors that fail at similar times. After visualizing that relationship, the template removes constant columns and linearly interpolates numeric measurements along the time index.

```r


# Summary of missing values
missing_summary <- colSums(is.na(data_df))
missing_summary <- data.frame(variable = names(missing_summary),
                              missing_count = missing_summary) %>%
  filter(missing_count > 0)

# Create a logical matrix to identify when values are missing
missing_matrix <- is.na(data_df)

# Investigate the correlation of missingness across sensors
missing_pattern <- as.data.frame(missing_matrix)

missing_pattern <- missing_pattern %>%
  select(-DateTime)

# Compute correlation missing_pattern 
missing_correlation <- cor(missing_pattern, use = "pairwise.complete.obs")
  

# Convert to long format for plotting
missing_correlation_melt <- melt(missing_correlation)
  
# Plot the correlation matrix
ggplot(missing_correlation_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                         name="Correlation") +
  theme_minimal() +
  labs(title = "Correlation of Missingness Across Sensors", 
       x = "Sensors",
       y = "Sensors") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  



# Convert tsibble to a data frame for easier manipulation
data_cleaned <- data_tsibble %>%
  as.data.frame()


# Identify and remove constant columns or non - numeric columns
constant_columns <- sapply(data_cleaned, function(x) length(unique(x[!is.na(x)])) == 1)

# Remove constant columns
data_cleaned <- data_cleaned[, !constant_columns ]

# Apply linear interpolation only on numeric columns
data_imputed <- data_cleaned %>%
  mutate(across(where(is.numeric), ~ na.interpolation(.)))


# Remove any rows with NA values after interpolation
data_imputed <- na.omit(data_imputed)

# Convert back to tsibble
data_imputed_tsibble <- data_imputed %>%
  as_tsibble(index = DateTime)


# Plotting the cleaned data
# Create a long format for plotting
data_long <- pivot_longer(data_imputed_tsibble, 
                          cols = -DateTime, 
                          names_to = "Variable", 
                          values_to = "Value")


# Create the plot
ggplot(data_long, aes(x = DateTime, y = Value, color = Variable)) +
  geom_line() +  
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +  # Create a separate panel for each variable
  labs(title = "Interpolation Air Quality Measurements Over Time",
       x = "DateTime",
       y = "Value") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




```

## Repeat PCA after interpolation

The imputed series support a second standardized PCA using more of the observed time span. The final steps inspect explained variance and scores, reconstruct the measurements from the first components, and compare reconstructed values with the cleaned series.

```r
  

# Convert to data frame for PCA
data_imputed_pca <- data_imputed_tsibble %>%
  as.data.frame() %>%
  select(-DateTime)


# Perform PCA
pca_result <- prcomp(data_imputed_pca, center = TRUE, scale. = TRUE)

# Check PCA results
summary(pca_result)


# Prepare data for the scree plot
scree_data <- data.frame(Variance = pca_result$sdev^2)  # Variance for each PC
scree_data$PC <- seq_along(scree_data$Variance)  # PC index


# Create Scree Plot
ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(1, nrow(scree_data)))  



  
# Biplot for PC1 vs PC2
biplot(pca_result, choices = c(1, 2), main = "PCA Biplot: PC1 vs PC2")

# Biplot for PC2 vs PC3
biplot(pca_result, choices = c(2, 3), main = "PCA Biplot: PC2 vs PC3")

# Biplot for PC3 vs PC4
biplot(pca_result, choices = c(3, 4), main = "PCA Biplot: PC3 vs PC4")

  
  
  
# Calculate the variance explained by each PC
explained_variance <- pca_result$sdev^2  # Variance for each PC
proportion_variance <- explained_variance / sum(explained_variance)  # Proportion of variance

# Cumulative variance explained
cumulative_variance <- cumsum(proportion_variance)

# Create a data frame for easier viewing
variance_data <- data.frame(
  PC = 1:length(cumulative_variance),
  Proportion = proportion_variance,
  Cumulative = cumulative_variance
)

# Loadings (rotation matrix)
loadings <- as.data.frame(pca_result$rotation)

# Choose the number of principal components (PCs) to retain for reconstruction
k <- min(7, ncol(pca_result$x))

# Project the original data into the PCA space (use only the first 'k' components)
pca_scores <- pca_result$x[, 1:k]

# Reconstruct the original data using the first 'k' PCs
reconstructed_scaled <- pca_scores %*% t(pca_result$rotation[, 1:k])
reconstructed_data <- sweep(reconstructed_scaled, 2, pca_result$scale, `*`)
reconstructed_data <- sweep(reconstructed_data, 2, pca_result$center, `+`)

# Ensure that both 'data_imputed_pca' and 'reconstructed_data' are matrices
# and have the same dimensions
# Calculate reconstruction error
reconstruction_error <- data_imputed_pca - reconstructed_data


# Combine the reconstructed data with original variables
final_data <- as.data.frame(reconstructed_data)
colnames(final_data) <- colnames(data_imputed_pca)

# Add back DateTime in the final dataset
final_data <- cbind(DateTime = data_imputed_tsibble$DateTime, final_data)


# Reshape for plotting
cleaned_data_long <- data_imputed_tsibble %>%
  as.data.frame() %>%
  pivot_longer(cols = -DateTime, names_to = "Variable", values_to = "Value")

reconstructed_data_long <- final_data %>%
  pivot_longer(cols = -DateTime, names_to = "Variable", values_to = "Value")


# Combine datasets for comparison
comparison_data <- bind_rows(
  cleaned_data_long %>% mutate(Type = "Cleaned"),
  reconstructed_data_long %>% mutate(Type = "Reconstructed")
)


# Extract PCA scores 
pca_scores <- as.data.frame(pca_result$x)

# PCA scores plot for PC1 vs PC2
ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.3) +  
  labs(title = "PCA Scores: PC1 vs PC2",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()


# Create a plot of the PCA scores for PC2 vs PC3
ggplot(pca_scores, aes(x = PC2, y = PC3)) +
  geom_point(alpha = 0.3) +  
  labs(title = "PCA Scores: PC2 vs PC3",
       x = "Principal Component 2",
       y = "Principal Component 3") +
  theme_minimal()


# PCA scores plot for PC3 vs PC4
ggplot(pca_scores, aes(x = PC3, y = PC4)) +
  geom_point(alpha = 0.3) + 
  labs(title = "PCA Scores: PC3 vs PC4",
       x = "Principal Component 3",
       y = "Principal Component 4") +
  theme_minimal()


# Plot individual plots for each variable
ggplot(comparison_data, aes(x = DateTime, y = Value, color = Type)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +  # Create a facet for each variable
  labs(title = "Comparison of Cleaned and Reconstructed Data by Variable",
       x = "DateTime",
       y = "Value") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Zooming example
# Combine PCA scores with DateTime for filtering
pca_scores_with_time <- cbind(DateTime = data_imputed_tsibble$DateTime, pca_scores)

# Filter for a specific time interval
filtered_scores <- pca_scores_with_time %>%
  filter(DateTime >= "2004-12-31" & DateTime < "2005-01-01")

# Plot the filtered PCA scores
ggplot(filtered_scores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5) +   # Add points
  geom_line(alpha = 0.5, color = "blue") +   # Add cumulative line connecting points
  labs(title = "PCA Scores (2004-12-31 to 2005-01-01): PC1 vs PC2",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()



  
  
  
  
```
