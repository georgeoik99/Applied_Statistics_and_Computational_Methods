# Computational Statistics

An R-first collection of MSc coursework on numerical optimization, simulation, Monte Carlo inference, and resampling. The notebooks retain the main examples and coding style of the original course material, with targeted corrections for reproducibility, parameterization, numerical stability, and interpretation.

## Notebooks

| Notebook | Main topics |
|---|---|
| [Numerical optimization](01_numerical_optimization.qmd) | Newton–Raphson, stopping rules, `optimize()`, `optim()`, Rosenbrock function |
| [Numerical MLE and likelihood inference](02_numerical_mle_and_likelihood_inference.qmd) | Numerical MLE, Hessian-based standard errors, Wald and likelihood-ratio intervals, full profile likelihood for parameters and functions of parameters |
| [Random generation and simulation](03_random_generation_and_simulation.qmd) | R distribution functions, inverse transform, rejection sampling, mixture and hierarchical simulation |
| [Monte Carlo methods](04_monte_carlo_methods.qmd) | Monte Carlo integration, estimator comparison, simulation-based hypothesis testing |
| [Bootstrap and resampling](05_bootstrap_and_resampling.qmd) | Non-parametric and parametric bootstrap, bootstrap-$t$, confidence-interval coverage, bias, standard errors, and dispersion checks |
| [Likelihood, bootstrap, and rejection sampling](06_likelihood_bootstrap_rejection_exercises.qmd) | Inverse-Gamma MLE, exponential-rate bootstrap, and Gamma rejection sampling exercises |
| [Gamma rate estimation](07_gamma_rate_estimation.qmd) | Newton–Raphson estimation and closed-form verification for a Gamma rate parameter |
| [Inverse-Gaussian inference and simulation](08_inverse_gaussian_inference_and_simulation.qmd) | Inverse-Gaussian likelihood, joint estimation, profile inference, bootstrap intervals, and random generation |
| [Monte Carlo variance reduction](09_monte_carlo_variance_reduction.qmd) | Antithetic variables, control variates, importance sampling, stratification, and empirical efficiency comparison |
| [Jackknife, BCa, and regression bootstrap](10_jackknife_bca_and_regression_bootstrap.qmd) | Jackknife bias and standard errors, BCa intervals, cases bootstrap, residual bootstrap, and influence diagnostics |
| [Permutation tests](11_permutation_tests.qmd) | Two-sample randomization tests, Monte Carlo p-values, energy distance, and distance correlation |
| [MCMC from scratch](12_mcmc_from_scratch.qmd) | Random-walk Metropolis, proposal tuning, effective sample size, multiple-chain diagnostics, and Gibbs sampling |

The `data/` directory contains the course datasets used by the notebooks.

## Reproduce the results

The statistical code uses base R and the datasets included in this folder. Rendering requires R, Quarto, `knitr`, and `rmarkdown`:

```r
install.packages(c("knitr", "rmarkdown"))
```

Then, from RStudio's Terminal:

```bash
quarto render
```

All stochastic examples set an explicit seed. Distribution parameters are labelled as `rate` or `scale` wherever the distinction matters.

### Work interactively in RStudio

1. Open `03_Computational_Statistics.Rproj`.
2. Open any `.qmd` notebook; it opens in the Visual editor by default.
3. Use the green **Run** button on a code chunk, or choose **Run All**.
4. Tables, numerical output, and plots appear directly below the corresponding R code.

The `.qmd` files are the editable notebooks. A full Render creates a local HTML version inside the ignored `_output/` directory; interactive chunk output remains visible directly below the code in RStudio.

## Academic context

The core examples originate from the MSc course **Computational Statistics** at the Athens University of Economics and Business. Course notes, classroom scripts, homework, and exercises were reviewed together; duplicated or methodologically unsafe fragments were consolidated rather than reproduced verbatim.

Notebooks 09–12 extend that core with selected topics from Maria L. Rizzo's *Statistical Computing with R* (2nd edition): variance reduction, jackknife and BCa inference, regression bootstrap, permutation methods, and introductory MCMC.
