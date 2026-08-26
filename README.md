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

Each analysis has a GitHub-native Markdown document and a matching standalone R script in the same directory. The `.md` file presents the explanation and source code directly on GitHub; the `.R` file follows the same order and includes the narrative as comments.

```text
CA2-Preprocessing/
  task1.md / task1.R        Air-quality preprocessing and PCA
  Task2.md / Task2.R        Frost temperature retrieval and diagnostics
CA3-Forecasting/
  Assignment03.md / Assignment03.R
CA4-Classification/
  Assignment_4.md / Assignment_4.R
```

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
Rscript CA3-Forecasting/Assignment03.R
```

The Frost template requires a client ID. Store it outside the repository and expose it only for the current shell:

```sh
export FROST_CLIENT_ID="your-client-id"
Rscript CA2-Preprocessing/Task2.R
```

Paths are resolved for execution from either the repository root or the relevant template directory. The analyses that use omitted datasets expect local copies beside their scripts.

## Data sources

- **UCI Air Quality:** Used by the preprocessing and PCA template and included in this repository. Source: UCI Machine Learning Repository, [Air Quality dataset](https://doi.org/10.24432/C59K5F), licensed under CC BY 4.0.
- **NMBU BIOKLIM weather data:** Daily temperature and irradiation used by the forecasting template. Local course copies are not redistributed; annual data are available from NMBU's [meteorological-data page](https://www.nmbu.no/forskning/grupper/meteorologiske-data). No redistribution licence is asserted here.
- **Elhub electricity consumption:** Aggregated municipal consumption used by the forecasting and classification templates. The local Ås subset is not redistributed; source data are available from Elhub's [data catalogue](https://elhub.no/data-og-innsikt/datakatalog) under CC BY 4.0.
- **ENTSO-E generation:** German generation by production type used by the classification template. The local export is not redistributed; use the ENTSO-E Transparency Platform's [Actual Generation per Production Type](https://transparency.entsoe.eu/generation/r2/actualGenerationPerProductionType/show?name=) view. The applicable open data are available under CC BY 4.0.
