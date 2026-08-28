# High-Dimensional Statistics

An R-first collection of MSc work on dimension reduction, multiple testing, classification, clustering, streaming regression, and count-data modeling. The notebooks retain the original straightforward coding style while applying targeted corrections for reproducibility and statistical interpretation.

## Notebooks

| Notebook | Main topics |
|---|---|
| [MNIST: PCA and t-SNE](01_mnist_pca_tsne.md) | Image inspection, PCA, t-SNE, common sampling, and careful interpretation of nonlinear embeddings |
| [Leukemia and multiple testing](02_leukemia_multiple_testing.md) | Gene-expression summaries, PCA, gene-wise t-tests, FWER, FDR, pFDR, and exploratory DE-gene visualization |
| [Streaming big-data regression](03_streaming_big_data_regression.md) | Chunk-wise cross-products, exact recursive OLS, random partitions, inverse-variance weighting, and estimator comparison |
| [Bike Sharing count models](04_bike_sharing_count_models.md) | Poisson, quasi-Poisson, negative binomial regression, overdispersion, interactions, and temporal evaluation |
| [Hotel booking classification](05_hotel_booking_classification.md) | Stratified evaluation, logistic regression, stepwise selection, Lasso, Random Forest, SVM, ROC AUC, and balanced accuracy |
| [Coffee clustering](06_coffee_clustering.md) | Hierarchical clustering, K-means, mclust, pgmm, fabMix, silhouette analysis, BIC, and post-hoc ARI |
| [Binary-data clustering](07_binary_data_clustering.md) | Asymmetric binary distance, missing values, Bernoulli mixtures, Bayesian mixtures, MDS, and cluster-number uncertainty |

The notebooks live directly in this directory; no additional `Assignment_1` subdirectory is required.

## Requirements

Core packages:

```r
install.packages(c(
  "ggplot2", "knitr", "MASS", "Rtsne", "keras3",
  "caret", "glmnet", "randomForest", "kernlab", "pROC",
  "cluster", "flexmix", "mclust"
))
keras3::install_keras(backend = "tensorflow")
```

Optional packages for the computationally intensive mixture analyses:

```r
install.packages(c("pgmm", "fabMix", "BayesBinMix"))
```

Bioconductor packages for the leukemia analysis:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("leukemiasEset", "qvalue"))
```

Install packages once from the R console, not inside the notebooks.

## Data

- MNIST is loaded through `keras3::dataset_mnist()`.
- The leukemia expression data are loaded from the Bioconductor `leukemiasEset` package.
- The regression dataset is generated reproducibly and processed through a temporary CSV without loading the whole file into memory.
- `data/hour.csv` is the Bike Sharing dataset from the [UCI Machine Learning Repository](https://archive.ics.uci.edu/dataset/275/bike%2Bsharing%2Bdataset), distributed under CC BY 4.0. Its original metadata are preserved in `data/bike_sharing_readme.txt`.
- `data/room_bookings.csv` is the 2,000-row course-provided sample used for the hotel cancellation assignment. The identifier and reservation date are excluded from modeling.
- `data/coffee.csv` reproduces the 43-sample `coffee` dataset distributed with the [CRAN `pgmm` package](https://CRAN.R-project.org/package=pgmm).
- `data/x1.txt` is the original 200-by-100 binary dataset supplied for the high-dimensional clustering exercise.

## Run in RStudio

Open any `.qmd` file and use **Run All**. Tables, numerical output, and plots will appear directly below their corresponding R chunks.

To perform a full render from the terminal:

```powershell
quarto render 01_mnist_pca_tsne.qmd
quarto render 02_leukemia_multiple_testing.qmd
quarto render 03_streaming_big_data_regression.qmd
quarto render 04_bike_sharing_count_models.qmd
quarto render 05_hotel_booking_classification.qmd
quarto render 06_coffee_clustering.qmd
quarto render 07_binary_data_clustering.qmd
```

The MNIST and leukemia notebooks require their external packages and data downloads before their first execution. The streaming-regression, Bike Sharing, hotel-classification, and core clustering analyses use local data. The `pgmm`, `fabMix`, and `BayesBinMix` sections are optional and controlled by flags near the top of their notebooks.

## Academic context

The material originates from the MSc course **High Dimensional Statistics** at the Athens University of Economics and Business. The original assignment logic and coding style are retained. Changes focus on correcting rejection counts, implementing exact recursive OLS, using genuine random partitions, removing response leakage, using train-only model selection, matching clustering methods with suitable distances, retaining all observations after explicit imputation, and treating MCMC boundary behavior cautiously. The social-network exercise from Assignment 2 is intentionally excluded.
