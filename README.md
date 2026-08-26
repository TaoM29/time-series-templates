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

Each analysis has a GitHub-native Markdown document and a matching standalone R script in the same directory. The `.md` file presents the explanation and source code directly on GitHub; the `.R` file preserves the same executable order with section comments for standalone use.

| Area | Documented analysis | Standalone script |
| --- | --- | --- |
| Preprocessing | [Air-quality preprocessing and PCA](preprocessing/air-quality-pca.md) | [`air-quality-pca.R`](preprocessing/air-quality-pca.R) |
| Preprocessing | [Daily temperature preprocessing with Frost](preprocessing/temperature-preprocessing.md) | [`temperature-preprocessing.R`](preprocessing/temperature-preprocessing.R) |
| Forecasting | [Electricity demand forecasting with weather covariates](forecasting/electricity-demand-forecasting.md) | [`electricity-demand-forecasting.R`](forecasting/electricity-demand-forecasting.R) |
| Clustering and classification | [Time-series clustering, segmentation, and classification](clustering-and-classification/time-series-clustering-and-classification.md) | [`time-series-clustering-and-classification.R`](clustering-and-classification/time-series-clustering-and-classification.R) |

The included UCI Air Quality data are stored beside the preprocessing templates as `preprocessing/AirQualityUCI.csv`. Other analyses resolve their locally downloaded source files from the relevant analysis directory.

The Markdown documents intentionally contain no saved plots, tables, console output, or execution state. Run the paired scripts locally to generate results. Generated reports and intermediate exports are ignored so the repository stays lightweight.

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

Read the `.md` files directly on GitHub, or open them in any Markdown viewer. Run a matching script from the repository root with, for example:

```sh
Rscript forecasting/electricity-demand-forecasting.R
```

The Frost template requires a client ID. Store it outside the repository and expose it only for the current shell:

```sh
export FROST_CLIENT_ID="your-client-id"
Rscript preprocessing/temperature-preprocessing.R
```

Paths are resolved for execution from either the repository root or the relevant template directory. The analyses that use omitted datasets expect local copies beside their scripts.

## Data sources

- **UCI Air Quality:** Used by the preprocessing and PCA template and included in this repository. Source: UCI Machine Learning Repository, [Air Quality dataset](https://doi.org/10.24432/C59K5F), licensed under CC BY 4.0.
- **NMBU BIOKLIM weather data:** Daily temperature and irradiation used by the forecasting template. Local archived copies are not redistributed; annual data are available from NMBU's [meteorological-data page](https://www.nmbu.no/forskning/grupper/meteorologiske-data). No redistribution licence is asserted here.
- **Elhub electricity consumption:** Aggregated municipal consumption used by the forecasting and clustering/classification templates. The local Ås subset is not redistributed; source data are available from Elhub's [data catalogue](https://elhub.no/data-og-innsikt/datakatalog) under CC BY 4.0.
- **ENTSO-E generation:** German generation by production type used by the clustering/classification template. The local export is not redistributed; use the ENTSO-E Transparency Platform's [Actual Generation per Production Type](https://transparency.entsoe.eu/generation/r2/actualGenerationPerProductionType/show?name=) view. The applicable open data are available under CC BY 4.0.
