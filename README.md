# Simulation and Data Analysis Codes for "Confidence Intervals and Hypothesis Testing for High-Dimensional Regression: A Review and Discussion"

This repository contains the replication codes for the report "Confidence Intervals and Hypothesis Testing for High-Dimensional Regression: A Review and Discussion" by Benjamin Smith. For the report, see [here](https://github.com/phd-statistics/sta3000-report/blob/main/report.pdf).

To ensure a smooth replication experience, please clone the repository and install the required packages listed below:

```
# Use the pak package to ensure installation is seemless
# install.packages("pak")
pak::pkg_install(
  c(
    "Matrix",
    "glmnet",
    "expm",
    "flare",
    "hdi",
    "mvnfast",
    "dplyr",
    "tidyr",
    "ggplot2",
    "patchwork" "doParallel",
    "foreach"
  )
)
```
# Simulation Study

To run the simulation study, ensure that you are in the in the `sims_and_da` directory and run the following in the R console.

```r
source("sim_main.R")
```

To post process the results from the analysis and create tables 1 and 3 in the paper, run the following in the R console:

```r
source("post.R")
```
# Data Analysis

To run the code which preforms the data analysis simply run in the R console.

```r
source("data_analysis.R")
```

# References

- [https://web.stanford.edu/~montanar/sslasso/](https://web.stanford.edu/~montanar/sslasso/)
