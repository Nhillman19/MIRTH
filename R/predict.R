#' Harmonize and generate predictions on out-of-sample data with MIRTH
#'
#' Applies a fitted [mirth()] object to new data whose outcomes are
#' unknown.
#' 1. The outcome and any missing covariates are multiply imputed with
#'    [mice::mice()], using the training data (where everything is observed)
#'    together with the new covariates and imaging features
#' 2. Each of the `M` completed data sets is harmonized with the batch-effect
#'    estimates learned on the *training* data, so no test information leaks
#'    into the harmonization model
#' 3. The prediction model is applied to each harmonized data set and the `M`
#'    predictions are pooled (mean for continuous outcomes; majority vote or
#'    averaged class probabilities for categorical outcomes).
#'
#' @details If `new_bat` contains batches that were not in the training data, their
#' batch effects must be estimated from the new data. MIRTH then alternates
#' between imputation and harmonization for up to `maxiter` rounds (stopping
#' early when the imputed values stabilize) before the final imputation. Note 
#' that this procedure is experimental and should be evaluated for multiple values
#' of `maxiter` if convergence does not occur
#'
#' @param object A fitted `mirth` object.
#' @param new_img Numeric matrix or data frame of imaging features for the new participants, with the same
#'   columns as the training data.
#' @param new_bat Vector of batch/site/scanner labels for the new participants.
#' @param new_covar Optional data frame of non-imaging covariates for the new participants, with the same columns
#'   as the training covariates. `NA` values are imputed. If `NULL` and the
#'   model used covariates, all covariates are imputed.
#' @param M Number of imputations.
#' @param imp_method Imputation method(s) passed to [mice::mice()] as `method`.
#'   A single string (e.g. `"rf"`, `"pmm"`,`"cart"`) applies to every incomplete variable; 
#'   a named character vector (names from the covariates and `"outcome"`) sets methods per variable.
#'   `NULL` uses the 'mice' defaults.
#' @param impute_with_img Logical variable indicating whether to use the imaging features as predictors in
#'   the imputation model. Defaults to TRUE.
#' @param mice_args Named list of further arguments to [mice::mice()], e.g.
#'   `list(maxit = 10)`.
#' @param maxiter Maximum number of imputation/harmonization iterations when unseen
#'   batches are present in the new data.
#' @param tol Convergence tolerance when iterating between imputation and harmonization. 
#' Compared against either the proportion of changed imputed classes, or the change 
#' in imputed continuous values relative to their standard deviation.
#' @param type For classification settings where `"response"` returns predicted values/classes and `"prob"` returns
#'   class probabilities averaged over imputations.
#' @param pool How categorical predictions are pooled across imputations, with options
#'   `"vote"` (majority vote) or `"prob"` (class with the highest average probability).
#' @param details Logical. If `TRUE`, returns a list with the pooled result
#'   together with per-imputation predictions, harmonized imaging data, and
#'   imputed covariates/outcomes.
#' @param seed Optional integer seed. The global random-number state is
#'   restored on exit.
#' @param ... Not used. Any arguments supplied here are ignored with a warning
#'
#' @return By default a numeric vector (regression), a factor
#'   (classification, `type = "response"`), or a matrix of class probabilities
#'   (`type = "prob"`). With `details = TRUE`, a list with elements `pred`,
#'   `prob` (classification), `pred_by_imputation`, `harmonized_img`,
#'   `imputed`, `unseen_batches`, and `iterations`.
#' @export

predict.mirth <- function(object, new_img, new_bat, new_covar = NULL, M = 10,
                          imp_method = "rf", impute_with_img = TRUE,
                          mice_args = list(), maxiter = 5, tol = 0.01,
                          type = c("response", "prob"), pool = c("vote", "prob"),
                          details = FALSE, seed = NULL,...) {
  if (...length()) {
    extra <- names(list(...))
    if (is.null(extra)) extra <- rep("", ...length())
    extra[extra == ""] <- "<unnamed>"
    warning("Ignoring unused argument(s) in `...`: ", paste(extra, collapse = ", "),
            call. = FALSE)
  }
  if (!is.null(seed)) withr::local_seed(seed)
  type <- match.arg(type)
  pool <- match.arg(pool)
  M <- .check_count(M, "M")
  maxiter <- .check_count(maxiter, "maxiter", min = 0)
  if (!is.list(mice_args)) stop("`mice_args` must be a list.", call. = FALSE)
  classif <- object$outcome_type == "classification"
  if (type == "prob" && !classif) {
    stop("`type = \"prob\"` is only available for categorical outcomes.", call. = FALSE)
  }

  tr <- object$train
  img <- .align_img(new_img, colnames(tr$img))
  n <- nrow(img)
  bat <- .prepare_bat(new_bat, n, "new_bat")
  covar <- .align_covar(new_covar, object$covar_template, n)

  unseen <- setdiff(levels(bat), levels(tr$bat))
  if (length(unseen)) {
    tiny <- unseen[table(bat)[unseen] < 2]
    if (length(tiny)) {
      stop("Batch effects for new batch(es) ", paste(tiny, collapse = ", "),
           " cannot be estimated from fewer than 2 participants.", call. = FALSE)
    }
    message("New batch(es) not seen in training: ", paste(unseen, collapse = ", "),
            ". Iterating between imputation and harmonization (maxiter = ",
            maxiter, ").")
  }

  impute <- function(train_img, new_img) {
    .impute(object, covar, if (impute_with_img) train_img, if (impute_with_img) new_img,
            n_new = n, M = M, imp_method = imp_method, mice_args = mice_args)
  }
  harmonize <- function(imp) {
    .as_img_matrix(
      object$harmonizer$transform(model = object$harm_model, newdata = img,
                                  newbat = bat, newcovar = imp),
      colnames(img), n, "transform"
    )
  }

  # Impute with original imaging data when every batch was seen in training
  imputed <- impute(tr$img, img)
  iterations <- 0L
  if (length(unseen) && maxiter > 0) {
    # New batches: refine the imputations in the harmonized feature space.
    for (it in seq_len(maxiter)) {
      current <- .collapse_imputations(imputed)
      img_cur <- harmonize(current)
      imputed <- impute(tr$img_harmonized, img_cur)
      iterations <- it
      if (.converged(current, .collapse_imputations(imputed), tol)) break
    }
  }

  # Harmonize the original new data once per imputation, then predict.
  harmonized <- lapply(imputed, harmonize)
  per_imp <- Map(function(h, imp) {
    use_cov <- object$model_covariates && !is.null(object$covar_template)
    nd <- .model_frame(h, if (use_cov) imp[names(object$covar_template)], object$feat_names)
    if (classif) {
      list(class = as.character(stats::predict(object$model, newdata = nd, type = "raw")),
           prob = as.matrix(stats::predict(object$model, newdata = nd, type = "prob")))
    } else {
      list(value = as.numeric(stats::predict(object$model, newdata = nd)))
    }
  }, harmonized, imputed)

  if (classif) {
    safe <- names(object$class_map)
    prob <- Reduce(`+`, lapply(per_imp, `[[`, "prob")) / M
    prob <- prob[, safe, drop = FALSE]
    votes <- do.call(cbind, lapply(per_imp, `[[`, "class"))
    pred_safe <- if (pool == "prob") safe[max.col(prob, ties.method = "first")]
                 else .vote(votes, prob, safe)
    pred <- factor(unname(object$class_map[pred_safe]), levels = object$outcome_levels)
    colnames(prob) <- object$outcome_levels
    by_imp <- matrix(object$class_map[votes], nrow = n)
    result <- if (type == "prob") prob else pred
  } else {
    by_imp <- do.call(cbind, lapply(per_imp, `[[`, "value"))
    pred <- rowMeans(by_imp)
    prob <- NULL
    result <- pred
  }

  if (!details) return(result)
  list(pred = pred, prob = prob, pred_by_imputation = by_imp,
       harmonized_img = harmonized, imputed = imputed,
       unseen_batches = unseen, iterations = iterations)
}

# ---- imputation ------------------------------------------------------------

# Returns a list of M complete data frames (one per imputation) for the new participants
.impute <- function(object, new_covar, train_img, new_img, n_new, M, imp_method,
                    mice_args) {
  tr <- object$train
  n_tr <- length(tr$outcome)
  new_outcome <- if (object$outcome_type == "classification") {
    factor(rep(NA, n_new), levels = object$outcome_levels)
  } else rep(NA_real_, n_new)

  dat <- data.frame(outcome = c(tr$outcome, new_outcome))
  if (!is.null(tr$covar)) dat <- cbind(rbind(tr$covar, new_covar), dat)
  if (!is.null(train_img)) {
    im <- as.data.frame(rbind(train_img, new_img))
    names(im) <- object$feat_names[colnames(train_img)]
    dat <- cbind(dat, im)
  }
  if (ncol(dat) < 2) {
    stop("Nothing to impute the outcome from: use `impute_with_img = TRUE` or ",
         "train with covariates.", call. = FALSE)
  }
  target <- setdiff(names(dat), object$feat_names)   # covariates + outcome
  incomplete <- names(dat)[vapply(dat, anyNA, logical(1))]

  method <- mice::make.method(dat)
  if (!is.null(imp_method)) {
    if (is.null(names(imp_method))) {
      if (length(imp_method) != 1) {
        stop("An unnamed `imp_method` must be a single string.", call. = FALSE)
      }
      method[incomplete] <- imp_method
    } else {
      bad <- setdiff(names(imp_method), target)
      if (length(bad)) {
        stop("`imp_method` names not found among covariates/outcome: ",
             paste(bad, collapse = ", "), call. = FALSE)
      }
      method[intersect(names(imp_method), incomplete)] <-
        imp_method[intersect(names(imp_method), incomplete)]
    }
  }
  args <- utils::modifyList(list(data = dat, m = M, method = method,
                                 printFlag = FALSE), mice_args)
  mids <- withCallingHandlers(
    do.call(mice::mice, args),
    warning = function(w) {
      if (grepl("logged events", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )

  new_rows <- n_tr + seq_len(n_new)
  keep <- c(names(object$covar_template), "outcome")
  out <- lapply(seq_len(M), function(m) {
    d <- mice::complete(mids, m)[new_rows, keep, drop = FALSE]
    rownames(d) <- NULL
    d
  })
  still_missing <- keep[vapply(keep, function(v) anyNA(out[[1]][[v]]), logical(1))]
  if (length(still_missing)) {
    stop("mice did not impute: ", paste(still_missing, collapse = ", "),
         ". This usually means mice dropped the variable as constant or collinear ",
         "(see mids$loggedEvents); try another `imp_method` or fewer predictors.",
         call. = FALSE)
  }
  out
}

# Collapse M imputations to one data set: mean for numeric, mode for factors.
.collapse_imputations <- function(imputed) {
  out <- imputed[[1]]
  for (v in names(out)) {
    vals <- lapply(imputed, `[[`, v)
    if (is.factor(out[[v]])) {
      mat <- do.call(cbind, lapply(vals, as.character))
      out[[v]] <- factor(apply(mat, 1, .mode), levels = levels(out[[v]]))
    } else {
      out[[v]] <- rowMeans(do.call(cbind, vals))
    }
  }
  out
}

.converged <- function(old, new, tol) {
  changes <- vapply(names(old), function(v) {
    a <- old[[v]]; b <- new[[v]]
    if (is.factor(a)) return(mean(as.character(a) != as.character(b)))
    s <- stats::sd(a)
    if (!is.finite(s) || s == 0) s <- 1
    mean(abs(a - b)) / s
  }, numeric(1))
  all(changes <= tol)
}

.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Majority vote across imputations; ties broken by the average probability.
.vote <- function(votes, prob, levels) {
  counts <- vapply(levels, function(l) rowSums(votes == l), numeric(nrow(votes)))
  counts <- matrix(counts, nrow = nrow(votes), dimnames = list(NULL, levels))
  counts <- counts + prob[, levels, drop = FALSE] * 1e-6
  levels[max.col(counts, ties.method = "first")]
}
