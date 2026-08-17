# Estimation and Statistical Inference

This directory contains Python-first notebooks that connect statistical theory with reproducible numerical examples, simulation, and interpretation. The material is organized as a compact learning sequence rather than a collection of disconnected code demonstrations.

The main topics are:

- likelihood and maximum likelihood estimation;
- finite-sample properties of estimators;
- ordinary least squares and Gaussian maximum likelihood;
- hypothesis testing and statistical power;
- sampling distributions and the Central Limit Theorem;
- probability distributions and Monte Carlo simulation;
- parametric and nonparametric testing.

## Contents

| Notebook | Main content | Status |
|---|---|---|
| `01_parameter_estimation.ipynb` | Likelihood, log-likelihood, MLE, bias, variance, standard error, MSE, exact confidence intervals, and simulation using the exponential-rate model | Implemented |
| `02_linear_regression_from_scratch.ipynb` | OLS derivation, Gaussian MLE, $\widehat{\beta}$, fitted values, residuals, RSS, error variance, coefficient covariance, standard errors, t-tests, F-test, confidence and prediction intervals, and diagnostics | Implemented |
| `03_hypothesis_testing.ipynb` | Null and alternative hypotheses, significance, p-values, Type I/II errors, power, parametric tests, and SciPy verification | Structured skeleton |
| `04_sampling_distributions_clt.ipynb` | Sampling distributions, CLT, normal approximation, Exercise 7 with a sum of Uniform variables, and Monte Carlo verification | Structured skeleton |
| `05_probability_distributions_and_simulation.ipynb` | PDF/PMF, CDF, quantiles, random generation, theoretical and empirical moments, and Monte Carlo examples | Implemented |
| `06_nonparametric_testing.ipynb` | Shapiro-Wilk, Mann-Whitney U, Wilcoxon signed-rank, Kruskal-Wallis, and parametric versus nonparametric interpretation | Implemented |
| `exercise_exponential_mle.ipynb` | Original standalone, executed version of Exercise 6 | Companion source notebook |

## Exercise 6 integration

Exercise 6 is **already integrated into `01_parameter_estimation.ipynb`** and serves as its central worked example. It assumes an independent sample from the exponential distribution parameterized by the rate $\theta$:

$$
f(x;\theta)=\theta e^{-\theta x}, \qquad x>0,\;\theta>0.
$$

The notebook derives

$$
\widehat{\theta}_{\mathrm{MLE}}
=\frac{n}{\sum_{i=1}^{n}X_i}
=\frac{1}{\bar X},
$$

then verifies the result numerically with `scipy.optimize` and studies its finite-sample bias, variance, standard error, MSE, confidence interval, and sampling distribution.

The separate `exercise_exponential_mle.ipynb` is intentionally retained as the original focused and executed exercise. It is not a second unresolved task and does not need to be repeated elsewhere.

## Suggested learning order

```text
01 Parameter Estimation
        ↓
02 Linear Regression from Scratch
        ↓
03 Hypothesis Testing
        ↓
04 Sampling Distributions and CLT
        ↓
05 Probability Distributions and Simulation
        ↓
06 Nonparametric Testing
```

The sequence is flexible. Notebook 05 can also be used as a probability and SciPy reference before beginning the inference notebooks.

## Environment

The implemented notebooks use:

- Python 3;
- NumPy;
- pandas;
- SciPy;
- Matplotlib.

Install the core dependencies with:

```bash
python -m pip install numpy pandas scipy matplotlib jupyter
```

Then open the directory in JupyterLab or VS Code and run the notebooks from top to bottom:

```bash
jupyter lab
```

All random experiments use explicit random-number generators or seeds so that results are reproducible.

## Scope and conventions

- Analytical results are introduced before library verification.
- Simulation supports the theory; it does not replace a mathematical argument.
- Parameterizations are stated explicitly, especially for exponential and gamma distributions.
- p-values are interpreted as evidence under the null model, not as the probability that the null hypothesis is true.
- Confidence intervals, effect estimates, assumptions, and practical interpretation are emphasized alongside significance tests.
- Regression metrics in this directory are educational, from-scratch calculations. Broader applied modeling workflows belong in `02_Statistical_Modeling/`.

## Academic context

The directory is based on MSc-level material in estimation, inference, probability distributions, hypothesis testing, and statistical computing. Examination exercises are expanded into reproducible Python studies while retaining the underlying mathematical derivations.

The notebooks are intended for learning and portfolio presentation. They should not be treated as a substitute for domain-specific statistical advice.
