#' Create your own harmonization method for MIRTH
#'
#' MIRTH allows you to use a custom harmonization model instead of the default
#' based on ComBat. For MIRTH, a harmonizer is a pair of functions -- one
#' that fits the harmonization model on training data and one that applies the
#' fitted model to new data. [harmonizer_comfam()] covers ComBat and its
#' extensions from 'ComBatFamily'; use `mirth_harmonizer()` to plug in any
#' other method that supports out-of-sample harmonization.
#'
#' @param fit A function with arguments `data` (numeric matrix, subjects x
#'   features), `bat` (factor), `covar` (data frame of covariates, including a
#'   column named `outcome`), and `formula` (a formula of the form
#'   `y ~ ...`). It must return a list with elements `model` (anything needed
#'   to harmonize new data) and `harmonized` (the harmonized training matrix,
#'   same dimensions as `data`).
#' @param transform A function with arguments `model` (the `model` element
#'   returned by `fit`), `newdata`, `newbat`, and `newcovar`, returning the
#'   harmonized version of `newdata` as a matrix. It must handle batches that
#'   were not seen during training, or raise an informative error.
#' @param name A short label used when printing.
#'
#' @return An object of class `mirth_harmonizer`.
#' @seealso [harmonizer_comfam()]
#' @export
#' @examples
#' # A "harmonizer" that does nothing, e.g. to benchmark against no harmonization
#' identity_harmonizer <- mirth_harmonizer(
#'   fit = function(data, bat, covar, formula) list(model = NULL, harmonized = data),
#'   transform = function(model, newdata, newbat, newcovar) newdata,
#'   name = "none"
#' )
mirth_harmonizer <- function(fit, transform, name = "custom") {
  if (!is.function(fit) || !is.function(transform)) {
    stop("`fit` and `transform` must both be functions.", call. = FALSE)
  }
  needed_fit <- c("data", "bat", "covar", "formula")
  needed_tr <- c("model", "newdata", "newbat", "newcovar")
  miss_fit <- setdiff(needed_fit, names(formals(fit)))
  miss_tr <- setdiff(needed_tr, names(formals(transform)))
  if ((length(miss_fit) && !"..." %in% names(formals(fit))) ||
      (length(miss_tr) && !"..." %in% names(formals(transform)))) {
    stop("`fit` must accept arguments ", paste(needed_fit, collapse = ", "),
         " and `transform` must accept ", paste(needed_tr, collapse = ", "), ".",
         call. = FALSE)
  }
  structure(list(name = name, fit = fit, transform = transform),
            class = "mirth_harmonizer")
}

#' ComBat-family harmonizer
#'
#' Harmonizes imaging features with [ComBatFamily::comfam()] and applies the
#' training-set estimates to new data with its `predict()` method. With the
#' default `model = stats::lm` this is standard ComBat; `model = mgcv::gam`
#' gives ComBat-GAM (use smooth terms such as `s(age)` in `base_harm_mod`).
#'
#' @param model Model function used to estimate covariate effects, passed to
#'   [ComBatFamily::comfam()].
#' @param eb Logical; use empirical Bayes shrinkage of site effects.
#' @param robust.LS Logical; use robust location and scale estimators.
#' @param ref_batch Optional reference batch; data from this batch (in both the
#'   training and new data) are left unchanged and other batches are mapped
#'   onto it.
#' @param ... Further arguments passed to `model`.
#'
#' @details CovBat ([ComBatFamily::covfam()]) is not offered because
#'   'ComBatFamily' does not currently provide out-of-sample harmonization for CovBat.
#'
#' @return An object of class `mirth_harmonizer`.
#' @export
#' @examples
#' harmonizer_comfam()
#' harmonizer_comfam(ref_batch = "site1")
harmonizer_comfam <- function(model = stats::lm, eb = TRUE, robust.LS = FALSE,
                              ref_batch = NULL, ...) {
  model <- match.fun(model)
  model_args <- list(...)
  force(eb)
  force(robust.LS)
  force(ref_batch)

  fit <- function(data, bat, covar, formula) {
    if (!is.null(ref_batch) && !ref_batch %in% levels(bat)) {
      stop("`ref_batch` (\"", ref_batch, "\") is not one of the training batches.",
           call. = FALSE)
    }
    obj <- do.call(ComBatFamily::comfam,
                   c(list(data = data, bat = bat, covar = covar, model = model,
                          formula = formula, eb = eb, robust.LS = robust.LS,
                          ref.batch = ref_batch),
                     model_args))
    list(model = obj, harmonized = obj$dat.combat)
  }
  transform <- function(model, newdata, newbat, newcovar) {
    out <- stats::predict(model, newdata = newdata, newbat = newbat,
                          newcovar = newcovar, robust.LS = robust.LS,
                          eb = eb)$dat.combat
    if (!is.null(ref_batch)) {
      ref <- as.character(newbat) == ref_batch
      if (any(ref)) out[ref, ] <- as.matrix(newdata)[ref, , drop = FALSE]
    }
    out
  }
  label <- if (identical(model, stats::lm)) "ComBat" else "ComBat family (custom model)"
  if (!is.null(ref_batch)) label <- paste0(label, ", reference batch = ", ref_batch)
  mirth_harmonizer(fit, transform, name = label)
}

#' @export
print.mirth_harmonizer <- function(x, ...) {
  cat("<mirth_harmonizer>", x$name, "\n")
  invisible(x)
}
