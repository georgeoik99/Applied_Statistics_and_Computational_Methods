# Zombie Apocalypse Logistic Regression

An educational statistical-modeling case study in R that uses fictional survey data to estimate the probability that an observation is classified as a zombie.

The project develops a complete, interpretable workflow:

- data validation and feature engineering;
- exploratory comparisons of humans and zombies;
- chi-square and Welch two-sample tests;
- binomial logistic regression with a logit link;
- coefficient, odds-ratio, confidence-interval, and p-value interpretation;
- likelihood-ratio testing and model diagnostics;
- in-sample classification metrics;
- probability predictions for fictional scenarios.

## Project structure

```text
01_Zombie_Apocalypse_Logistic_Regression/
├── README.md
├── zombie_apocalypse_logistic_regression.qmd
└── data/
    └── zombies.csv
```

After rendering for GitHub, Quarto also creates the corresponding Markdown output and figure assets.

## Statistical methods

The primary model is

\[
\operatorname{logit}\{P(Y=\text{Zombie}\mid X)\}
= \beta_0 + X^T\beta,
\]

fitted with R's `glm(..., family = binomial(link = "logit"))`.

The analysis reports:

- exploratory density and proportion plots;
- bivariate association tests;
- coefficient estimates and Wald tests;
- odds ratios with 95% confidence intervals;
- an overall likelihood-ratio test against the intercept-only model;
- generalized variance-inflation factors;
- residual and linearity-in-the-logit diagnostics;
- accuracy, sensitivity, specificity, Brier score, and AIC.

Classification metrics are explicitly described as **in-sample** and should not be interpreted as estimates of out-of-sample performance.

## Requirements

- R 4.3 or later
- Quarto
- R packages: `ggplot2`, `gridExtra`, `car`, and `knitr`

Install the packages once from an R console:

```r
install.packages(c("ggplot2", "gridExtra", "car", "knitr"))
```

## Render

Open a terminal in this directory and run:

```powershell
quarto render zombie_apocalypse_logistic_regression.qmd
```

The document uses the GitHub-Flavored Markdown format so that code, results, tables, equations, and plots can be displayed directly on GitHub.

## Academic context

This project began as a guided DataCamp educational case study and was subsequently reorganized and expanded with clearer statistical explanations, diagnostics, reproducible paths, odds-ratio inference, model evaluation, and personal scenario predictions.

The dataset contains 200 fictional observations. The results are educational and must not be interpreted as real epidemiological or survival-risk estimates.

## Limitations

- The data and outcome are fictional.
- Bivariate tests are exploratory and do not replace multivariable reasoning.
- Apparent classification performance is evaluated on the training observations.
- The project does not include external validation or causal identification.
- Dataset redistribution rights should be confirmed before public publication.
