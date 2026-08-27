# Applied Statistics & Computational Methods

A collection of MSc coursework, statistical case studies, programming exercises and applied projects developed using **R** and **Python**.

The repository focuses on the practical implementation of statistical theory, mathematical optimization and computational methods. It includes academic assignments, classroom exercises and selected extended case studies from the MSc in Applied Statistics and Data Analytics.

## Current Repository Map

- [`01_Estimation_and_Inference/`](01_Estimation_and_Inference/) — Python notebooks on estimation, inference, distributions, regression from scratch and hypothesis testing.
- [`02_Statistical_Modeling/`](02_Statistical_Modeling/) — R/Quarto projects on linear models, generalized linear models, model selection, prediction, and diagnostics.
- [`03_Computational_Statistics/`](03_Computational_Statistics/) — R/Quarto case studies on numerical optimization, MLE, simulation, Monte Carlo and bootstrap methods.
- [`04_High_Dimensional_Statistics/`](04_High_Dimensional_Statistics/) — R/Quarto notebooks on PCA, t-SNE, multiple testing, streaming regression and high-dimensional applications.

## Main Areas

### Statistical Modeling

- Simple and multiple linear regression
- Model specification and interpretation
- Variable selection
- Interaction and polynomial terms
- Residual analysis and model diagnostics
- Multicollinearity and regularization
- Ridge regression
- Prediction and model evaluation

### Generalized Linear Models

Applications of the Generalized Linear Model framework for different types of response variables:

- Gaussian regression
- Binomial and logistic regression
- Poisson regression
- Negative binomial regression
- Gamma regression
- Inverse Gaussian models
- Multinomial response models
- Ordinal response models
- Count-data models
- Overdispersion analysis
- Zero-inflated model extensions

Selected exercises may also use **Vector Generalized Linear Models (VGLMs)** through the `VGAM` framework in R.

### Statistical Inference and Estimation

- Likelihood and log-likelihood functions
- Maximum Likelihood Estimation
- Method of Moments
- Least Squares Estimation
- Properties of estimators
- Bias, variance and consistency
- Confidence intervals
- Standard errors
- Asymptotic inference

### Hypothesis Testing

- Null and alternative hypotheses
- p-values and significance levels
- One-sample and two-sample t-tests
- Paired t-tests
- Proportion tests
- Chi-square tests
- Analysis of Variance
- Likelihood-ratio tests
- Parametric and non-parametric testing
- Multiple-testing considerations

### Computational Statistics

- Monte Carlo simulation
- Random-variable generation
- Numerical experiments
- Bootstrap estimation
- Bootstrap confidence intervals
- Parametric and non-parametric bootstrap
- Resampling methods
- Permutation-based procedures
- Numerical approximation of statistical quantities
- Simulation-based evaluation of estimators

### Optimization Methods

- Objective functions, gradients and Hessian matrices
- Gradient Descent
- Newton methods
- Quasi-Newton methods
- Symmetric Rank-One updates
- AdaGrad and adaptive optimization
- Trust-region methods
- Projected Gradient Descent
- Constrained optimization
- Lagrange multipliers and KKT conditions
- Numerical optimization for Ridge Regression

### Additional Topics

The repository may gradually include material related to:

- Multivariate statistical analysis
- Time-series analysis
- Bayesian statistics
- Big Data methods
- Statistical learning
- Sampling methods
- Experimental design
- Data visualization
- Statistical computing and algorithm implementation

## Repository Content

Projects will be organized by statistical area or university course. Depending on the case study, each folder may contain:

- R scripts or R Markdown files
- Python scripts or Jupyter notebooks
- Mathematical derivations
- Data preparation and exploratory analysis
- Statistical model implementation
- Diagnostic plots and results
- Assignment reports
- Supporting datasets or data-source instructions
- Project-specific documentation

Example structure:

```text
Applied_Statistics_and_Computational_Methods/
│
├── 01_Estimation_and_Inference/
├── 02_Statistical_Modeling/
├── 03_Computational_Statistics/
├── 04_High_Dimensional_Statistics/
└── README.md
```

The final structure will evolve as additional coursework and case studies are reviewed and added.

## Tools and Technologies

- R
- RStudio
- Python
- Jupyter Notebook
- NumPy
- Pandas
- SciPy
- Matplotlib
- Seaborn
- Base R
- tidyverse
- ggplot2
- VGAM

The exact dependencies will be documented separately within each project.

## Academic Context

A significant part of this repository originates from coursework completed during the **MSc in Applied Statistics and Data Analytics at the Athens University of Economics and Business**.

Projects will be clearly identified as:

- MSc coursework
- Classroom exercises
- Extended academic case studies
- Personal implementations or portfolio extensions

Course-provided code or methodological templates will be acknowledged where applicable. Any personal extensions, dataset selection, experimentation, analysis and interpretation will also be described explicitly.
