# time-series-templates

Reusable time-series analysis templates primarily written in R. The material began as DAT320 coursework and study work and has been cleaned, corrected, and organized as reference implementations rather than a course report.

## Topics

- Time-indexed data import, time-zone handling, gap detection, and regularization
- Missing-value analysis and linear interpolation
- Air-quality sensor exploration, principal component analysis (PCA), and reconstruction
- Stationarity checks with ADF and KPSS tests
- ACF/PACF analysis, seasonal differencing, and STL decomposition
- Rolling-origin cross-validation and ETS forecasting
- Granger predictability tests and ARIMAX models with weather covariates
- Time-series dissimilarity with Euclidean, correlation, AR.PIC, and dynamic time warping measures
- Hierarchical clustering, weekday/weekend classification, and hidden Markov model segmentation

## Organization

Each retained notebook has a matching standalone R script in the same directory. The scripts follow the notebook order and include the Markdown narrative as comments.

```text
CA2-Preprocessing/
  task1.Rmd / task1.R       Air-quality preprocessing and PCA
  Task2.Rmd / Task2.R       Frost temperature retrieval and diagnostics
CA3-Forecasting/
  Assignment03.Rmd / Assignment03.R
CA4-Classification/
  Assignment_4.Rmd / Assignment_4.R
```

The notebooks intentionally contain no saved plots, tables, console output, or execution state. Run them locally to generate results. Generated reports and intermediate exports are ignored so the repository stays lightweight.

## Requirements and use

Use a current R installation and install the packages needed by the template you plan to run. Across the repository these include:

```r
install.packages(c(
  "caret", "depmixS4", "dplyr", "dtw", "factoextra",
  "forecast", "ggplot2", "imputeTS", "jsonlite", "lmtest",
  "lubridate", "pheatmap", "readr", "readxl", "reshape2", "stringr",
  "tidyr", "TSclust", "tsibble", "tseries", "urca", "viridis", "zoo"
))
```

Open an `.Rmd` file in RStudio or another R Markdown environment to work interactively. Run a matching script from the repository root with, for example:

```sh
Rscript CA3-Forecasting/Assignment03.R
```

The Frost template requires a client ID. Store it outside the repository and expose it only for the current shell:

```sh
export FROST_CLIENT_ID="your-client-id"
Rscript CA2-Preprocessing/Task2.R
```

Input datasets are kept beside the analyses that use them. Paths are resolved for execution from either the repository root or the relevant template directory.
