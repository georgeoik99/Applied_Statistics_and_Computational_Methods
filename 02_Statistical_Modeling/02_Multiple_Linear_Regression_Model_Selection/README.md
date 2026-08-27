# Multiple Linear Regression and Model Selection

An R/Quarto MSc project examining a regression problem with one response and 52 candidate predictors. The analysis preserves the structure of the original assignment while making the workflow reproducible and correcting the broken prediction section.

## Methods

- descriptive statistics and distribution plots;
- grouped correlation analysis;
- full and screened multiple linear regression models;
- variance inflation factors;
- stepwise selection using AIC and BIC;
- coefficient estimates and confidence intervals;
- residual, leverage, and Cook's-distance diagnostics;
- fitted values and prediction intervals.

## Project structure

```text
02_Multiple_Linear_Regression_Model_Selection/
├── README.md
├── multiple_linear_regression_model_selection.qmd
└── data/
    └── multiple_regression_dataset2.csv
```

The CSV contains the 59 complete observations and 53 modeling variables used by the original assignment: response `Y2` and predictors `X1`–`X52`.

## Requirements

- R 4.3 or later
- Quarto
- R packages: `ggplot2`, `tidyr`, `corrplot`, `broom`, `car`, and `knitr`

Install the packages once from an R console:

```r
install.packages(c("ggplot2", "tidyr", "corrplot", "broom", "car", "knitr"))
```

## Render

From this directory:

```powershell
quarto render multiple_linear_regression_model_selection.qmd
```

The notebook also supports interactive execution in RStudio, where code output and plots appear below each chunk.

## Methodological note

The full model has 52 predictors but only 59 observations. It is retained because it belongs to the original coursework, but its in-sample fit is not treated as evidence of predictive accuracy. The final discussion explicitly addresses overfitting, multicollinearity, stepwise-selection instability, and post-selection inference.

## Academic context

This project originated as MSc coursework in linear models and was reorganized into a self-contained portfolio analysis. The original exploratory, VIF, AIC/BIC, and diagnostic workflow remains recognizable; changes are limited to reproducibility, error correction, presentation, and statistically necessary caveats.
