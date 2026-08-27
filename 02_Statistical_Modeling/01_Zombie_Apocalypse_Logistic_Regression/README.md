# Zombie Apocalypse Logistic Regression

A playful R project that uses fictional data to estimate the probability that a person is classified as a zombie.

## Methods

- exploratory comparisons of humans and zombies;
- chi-square and Welch two-sample tests;
- logistic regression;
- odds ratios and 95% confidence intervals;
- likelihood-ratio testing;
- confusion matrix and simple classification metrics;
- VIF, residual, linearity, and influence diagnostics;
- probability predictions for fictional scenarios.

## Project structure

```text
01_Zombie_Apocalypse_Logistic_Regression/
├── README.md
├── zombie_apocalypse_logistic_regression.qmd
└── data/
    └── zombies.csv
```

## Requirements

- R
- Quarto
- `ggplot2`
- `gridExtra`
- `car`

Install the packages once:

```r
install.packages(c("ggplot2", "gridExtra", "car"))
```

Open `zombie_apocalypse_logistic_regression.qmd` in RStudio and use **Run All**, or render it with:

```powershell
quarto render zombie_apocalypse_logistic_regression.qmd
```

## Origin and limitations

The project began as a guided DataCamp exercise and was later revised for clearer explanations, reproducible paths, corrected tests, model diagnostics, and probability predictions.

The dataset contains 200 fictional observations. The classification results are in-sample and the analysis is not a real medical, epidemiological, or survival-risk model.
