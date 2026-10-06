# MIRTH

**M**ultiple **I**mputation for **R**emoving **T**echnical **H**eterogeneity proposes
a leakage-free integration of ComBat harmonization into machine-learning
pipelines. When imaging data is pooled across sites, imaging features often contain *site effects*. Harmonizing the
whole data set (with training and test outcomes in the model) leaks test information into
training. Additionally, harmonizing without the outcome in the model removes true biological signal
when the outcome is imbalanced across sites. MIRTH aims to address these issues by

1. Fitting the harmonization model (with the outcome) on the training data only
2. Training a prediction model on the harmonized training data, with
   hyperparameters tuned by inner cross-validation
3. Imputing the unknown test outcomes (and any other missing test covariates) via multiple imputation
4. Harmonizing each imputed test set with the training-set batch-effect
   estimates
5. Applying the fitted prediction model to get a prediction for each harmonized test set 
6. Pooling predictions across imputations to achieve a final prediction

## Installation

MIRTH depends on [ComBatFamily](https://github.com/andy1764/ComBatFamily),
which is installed automatically from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Nhillman19/MIRTH")
library(MIRTH)
```

## Quick start

```r
# Toy data: 3 sites with additive site effects
set.seed(1)
n <- 150
site <- factor(sample(paste0("site", 1:3), n, replace = TRUE))
y <- factor(rbinom(n, 1, c(0.2, 0.5, 0.8)[site]), labels = c("control", "case"))
img <- matrix(rnorm(n * 10), n, 10) + 0.25*as.integer(site) + 0.25*(y == "case")
covar <- data.frame(age = rnorm(n, 70, 5))
train <- 1:120
test <- 121:150

fit <- mirth(img[train, ], y[train], site[train], covar[train, , drop = FALSE],
             model_method = "glm", inner_folds = 3, seed = 1)
fit
pred <- predict(fit, img[test, ], site[test], covar[test, , drop = FALSE],
                M = 3, imp_method = "pmm", seed = 1)
mean(pred == y[test])
```

## Features

| | |
|---|---|
| Outcomes | continuous, binary, multi-class (any class labels) |
| Missing data | test outcomes always; test covariates (`NA`, or omit `new_covar` entirely) |
| Prediction models | anything in [caret](https://topepo.github.io/caret/available-models.html) (`"rf"`, `"glmnet"`, `"svmRadial"`, ...) |
| Imputation | any [mice](https://amices.org/mice/) method, globally or per variable (`imp_method = c(outcome = "rf", age = "pmm")`) |
| Harmonization | ComBat, ComBat-GAM (`harmonizer_comfam(model = mgcv::gam)`), reference-batch ComBat (`harmonizer_comfam(ref_batch = "siteA")`), or your own via `mirth_harmonizer()` |
| New sites | (Experimental) batches absent from training are handled by alternating imputation and harmonization (`maxiter`, `tol`) |
| Diagnostics | `predict(..., details = TRUE)` returns per-imputation predictions, harmonized test data, and imputed values |

## Using your own harmonization method

Any method that can harmonize new data with parameters learned on training
data can be plugged in:

```r
my_harmonizer <- mirth_harmonizer(
  fit = function(data, bat, covar, formula) {
    # ... fit on training data ...
    list(model = fitted_object, harmonized = harmonized_training_matrix)
  },
  transform = function(model, newdata, newbat, newcovar) {
    # ... apply `model` to new data; return a matrix ...
  },
  name = "my method"
)
fit <- mirth(..., harmonizer = my_harmonizer)
```

## Acknowledgements

If you use the MIRTH package, please cite the original manuscript where the method was proposed:

Hillman, N., Chen, A., Hu, F., Vandekar, S., Melhem, R., Beason-Held, L., Satterthwaite, T., Davatzikos, C., Shou, H., & Shinohara, R. (2026). 
Imputation-Based Harmonization Mitigates Site Effects Without Data Leakage in Machine Learning Studies. bioRxiv, 
2026.09.28.755204. https://doi.org/10.64898/2026.09.28.755204

Since publication, the code has been reorganized to increase its modularity and 
robustness to argument specification, with Claude Opus 5.5 (Anthropic) providing valuable feedback. 
The updated code generates identical harmonized data and predictions when compared to the original version.

