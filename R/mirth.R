#' MIRTH: A Harmonization Procedure for Prediction Pipelines
#'
#' Fits a harmonization model on the training imaging data, harmonizes
#' the training data, and fits a prediction model on the harmonized features
#' with `caret`. Use [predict.mirth()] to harmonize and predict on held-out target data. The
#' unknown outcomes (and any missing covariates) are multiply imputed, while each
#' imputed data set is harmonized with the *training* batch-effect estimates,
#' and predictions are pooled across imputations.
#'
#' @param train_img Numeric matrix or data frame of imaging features
#'   (participants x features). Data must be fully observed (i.e., no NA entries)
#' @param train_outcome Outcome to predict represented as a vector. Numeric for regression; factor,
#'   character, or logical for classification (binary or multi-class).
#' @param train_bat Vector of batch/site/scanner labels for each participant.
#' @param train_covar Optional data frame of non-imaging covariates (fully observed
#'   in the training data). Cannot include variable names `y`,
#'   `batch`, or `outcome`, since these are reserved by MIRTH
#' @param base_harm_mod Covariate model for harmonization, as a formula or
#'   string, e.g. `y ~ age + sex` or `~ age + sex`. The outcome is added
#'   automatically as `+ outcome`. To control how it enters (e.g. `s(outcome)`
#'   with a GAM harmonizer), refer to it as `outcome` yourself. `NULL` (the
#'   default) uses all columns of `train_covar` in a linear model
#' @param harmonizer A [mirth_harmonizer()], which defaults to [harmonizer_comfam()]
#'   (ComBat).
#' @param model_method A `caret` model name (see `caret::modelLookup()`),
#'   e.g. `"rf"`, `"glmnet"`, `"svmRadial"`, `"lm"`,`"glm"`.
#' @param tune_grid Optional data frame of tuning parameters for `caret`.
#' @param inner_folds Number of cross-validation folds used for tuning.
#' @param metric Performance metric used to select tuning parameters. Defaults
#'   to `"RMSE"` (regression), `"ROC"` (binary), or `"logLoss"` (multi-class).
#' @param model_covariates Logical, indicates whether `train_covar` should be included as 
#' predictors in the prediction model in addition to imaging features.
#' @param train_args Named list of further arguments passed to
#'   [caret::train()], e.g. `list(ntree = 500)` or `list(trControl = ...)` to
#'   override the default cross-validation scheme.
#' @param seed Optional integer seed. The global random-number state is
#'   restored on exit.
#'
#' @return An object of class `mirth` containing the fitted harmonizer, the
#'   `caret` model, and the training data needed to impute new data.
#'
#' @references Hillman, N., Chen, A., Hu, F., Vandekar, S., Melhem, R., Beason-Held, L., Satterthwaite, T., Davatzikos, C., Shou, H., & Shinohara, R. (2026). 
#' Imputation-Based Harmonization Mitigates Site Effects Without Data Leakage in Machine Learning Studies. bioRxiv, 
#' 2026.09.28.755204. https://doi.org/10.64898/2026.09.28.755204
#'
#' @seealso [predict.mirth()], [harmonizer_comfam()]
#' @export
#' @examples
#' # Toy data: 3 sites with different case prevalence and additive site effects
#' set.seed(1)
#' n <- 150
#' site <- factor(sample(paste0("site", 1:3), n, replace = TRUE))
#' y <- factor(rbinom(n, 1, c(0.2, 0.5, 0.8)[site]), labels = c("control", "case"))
#' img <- matrix(rnorm(n * 10), n, 10) + 0.25*as.integer(site) + 0.25*(y == "case")
#' covar <- data.frame(age = rnorm(n, 70, 5))
#' train <- 1:120
#' test <- 121:150
#'
#' fit <- mirth(img[train, ], y[train], site[train], covar[train, , drop = FALSE],
#'              model_method = "glm", inner_folds = 3, seed = 1)
#' fit
#' pred <- predict(fit, img[test, ], site[test], covar[test, , drop = FALSE],
#'                 M = 3, imp_method = "pmm", seed = 1)
#' mean(pred == y[test])
mirth <- function(train_img, train_outcome, train_bat, train_covar = NULL,
                  base_harm_mod = NULL, harmonizer = harmonizer_comfam(),
                  model_method = "rf", tune_grid = NULL, inner_folds = 5,
                  metric = NULL, model_covariates = TRUE, train_args = list(),
                  seed = NULL) {
  if (!is.null(seed)) withr::local_seed(seed)
  if (!inherits(harmonizer, "mirth_harmonizer")) {
    stop("`harmonizer` must be created with harmonizer_comfam() or mirth_harmonizer().",
         call. = FALSE)
  }
  inner_folds <- .check_count(inner_folds, "inner_folds", min = 2)
  if (!is.list(train_args)) stop("`train_args` must be a list.", call. = FALSE)

  img <- .prepare_img(train_img, "train_img")
  n <- nrow(img)
  covar <- .prepare_covar(train_covar, n, "train_covar")
  if (!is.null(covar) && anyNA(covar)) {
    stop("`train_covar` contains missing values; MIRTH requires complete ",
         "training covariates (impute or drop them before fitting).", call. = FALSE)
  }
  bat <- .prepare_bat(train_bat, n, "train_bat")
  small <- table(bat) < 2
  if (any(small)) {
    stop("Each training batch needs at least 2 participants; too small batches: ",
         paste(names(small)[small], collapse = ", "), call. = FALSE)
  }
  outcome <- .prepare_outcome(train_outcome, n)

  # Harmonization model fit on training data
  harm_formula <- .build_harm_formula(base_harm_mod, names(covar), parent.frame())
  harm_fit <- harmonizer$fit(data = img, bat = bat,
                             covar = .harm_covar(covar, outcome$value),
                             formula = harm_formula)
  img_harm <- .as_img_matrix(harm_fit$harmonized, colnames(img), n, "fit")

  # Prediction model tuned by cross-validation on harmonized data
  feat_names <- .feature_names(colnames(img), names(covar))
  model <- .fit_prediction_model(
    img = img_harm, covar = if (model_covariates) covar else NULL,
    outcome = outcome, feat_names = feat_names, model_method = model_method,
    tune_grid = tune_grid, inner_folds = inner_folds, metric = metric,
    train_args = train_args
  )

  structure(
    list(
      call = match.call(),
      outcome_type = outcome$type,
      outcome_levels = outcome$levels,
      harmonizer = harmonizer,
      harm_model = harm_fit$model,
      harm_formula = harm_formula,
      model = model$fit,
      class_map = model$class_map,
      feat_names = feat_names,
      model_covariates = model_covariates,
      covar_template = if (is.null(covar)) NULL else covar[0, , drop = FALSE],
      train = list(img = img, img_harmonized = img_harm, covar = covar,
                   outcome = outcome$value, bat = bat)
    ),
    class = "mirth"
  )
}

# ---- helpers -------------------------------------------------------------

.build_harm_formula <- function(base, covar_names, env = parent.frame()) {
  if (is.null(base)) {
    rhs <- c(covar_names, "outcome")
    f <- stats::reformulate(rhs, response = "y")
  } else {
    if (is.character(base)) {
      if (length(base) != 1) stop("`base_harm_mod` must be a single string.", call. = FALSE)
      base <- stats::as.formula(base, env = env)
    }
    if (!inherits(base, "formula")) {
      stop("`base_harm_mod` must be a formula, a string, or NULL.", call. = FALSE)
    }
    if (length(base) == 2) base <- stats::as.formula(call("~", quote(y), base[[2]]), env = env)
    if (!identical(base[[2]], quote(y))) {
      stop("The left-hand side of `base_harm_mod` must be `y` (e.g. y ~ age + sex).",
           call. = FALSE)
    }
    vars <- all.vars(base[[3]])
    unknown <- setdiff(vars, c(covar_names, "outcome"))
    if (length(unknown)) {
      stop("`base_harm_mod` refers to variable(s) not in `train_covar`: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    f <- if ("outcome" %in% vars) base else stats::update(base, . ~ . + outcome)
  }
  environment(f) <- env
  f
}

.harm_covar <- function(covar, outcome) {
  if (is.null(covar)) data.frame(outcome = outcome) else data.frame(covar, outcome = outcome)
}

.as_img_matrix <- function(x, names, n, where) {
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) != n || ncol(x) != length(names)) {
    stop("The harmonizer's `", where, "` step returned data with the wrong ",
         "dimensions.", call. = FALSE)
  }
  if (anyNA(x)) stop("The harmonizer returned missing values.", call. = FALSE)
  colnames(x) <- names
  rownames(x) <- NULL
  x
}

# Make collision-free and valid names for imaging features
.feature_names <- function(img_names, covar_names) {
  all <- make.names(c(covar_names, "outcome", img_names), unique = TRUE)
  stats::setNames(all[-seq_len(length(covar_names) + 1)], img_names)
}

.model_frame <- function(img, covar, feat_names) {
  x <- as.data.frame(img)
  names(x) <- feat_names[colnames(img)]
  if (is.null(covar)) x else cbind(covar, x)
}

.fit_prediction_model <- function(img, covar, outcome, feat_names, model_method,
                                  tune_grid, inner_folds, metric, train_args) {
  df <- .model_frame(img, covar, feat_names)
  class_map <- NULL
  if (outcome$type == "classification") {
    # Class labels need to be valid R names
    safe <- make.names(outcome$levels, unique = TRUE)
    class_map <- stats::setNames(outcome$levels, safe)
    df$outcome <- factor(safe[as.integer(outcome$value)], levels = safe)
    binary <- length(safe) == 2
    summary_fun <- if (binary) caret::twoClassSummary else caret::mnLogLoss
    if (is.null(metric)) metric <- if (binary) "ROC" else "logLoss"
    ctrl <- caret::trainControl(method = "cv", number = inner_folds,
                                classProbs = TRUE, summaryFunction = summary_fun)
  } else {
    df$outcome <- outcome$value
    if (is.null(metric)) metric <- "RMSE"
    ctrl <- caret::trainControl(method = "cv", number = inner_folds)
  }
  args <- utils::modifyList(
    list(data = df, method = model_method, trControl = ctrl, metric = metric),
    train_args
  )
  if (!is.null(tune_grid)) args$tuneGrid <- tune_grid
  # the formula must be the first (unnamed) argument for S3 dispatch
  fit <- do.call(caret::train, c(list(outcome ~ .), args))
  list(fit = fit, class_map = class_map)
}

#' @export
print.mirth <- function(x, ...) {
  tr <- x$train
  cat("MIRTH prediction pipeline\n")
  cat("  Outcome:       ", if (x$outcome_type == "regression") "continuous" else
    paste0("categorical (", paste(x$outcome_levels, collapse = ", "), ")"), "\n", sep = "")
  n_cov <- if (is.null(tr$covar)) 0 else ncol(tr$covar)
  plural <- function(k, word) paste0(k, " ", word, if (k != 1) "s")
  cat("  Training data: ", plural(nrow(tr$img), "participant"), ", ",
      plural(ncol(tr$img), "imaging feature"), ", ", plural(n_cov, "covariate"),
      ", ", plural(nlevels(tr$bat), "site"), "\n", sep = "")
  cat("  Harmonizer:    ", x$harmonizer$name, "\n", sep = "")
  cat("  Harm. formula: ", paste(deparse(x$harm_formula, width.cutoff = 500), collapse = " "),
      "\n", sep = "")
  best <- x$model$bestTune
  cat("  Model:         ", x$model$method,
      if (!is.null(best) && !identical(names(best), "parameter"))
        paste0(" (", paste(names(best), best, sep = " = ", collapse = ", "), ")"),
      "\n", sep = "")
  res <- x$model$results
  m <- x$model$metric
  if (!is.null(res) && m %in% names(res)) {
    val <- if (x$model$maximize) max(res[[m]], na.rm = TRUE) else min(res[[m]], na.rm = TRUE)
    cat("  Inner CV ", m, ":  ", format(val, digits = 3), "\n", sep = "")
  }
  invisible(x)
}
