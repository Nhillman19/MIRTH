# Input checking and standardization for each aspect of the provided data

.reserved_names <- c("y", "batch", "outcome")

.prepare_img <- function(img, arg) {
  if (is.data.frame(img)) {
    is_num <- vapply(img, is.numeric, logical(1))
    if (!all(is_num)) {
      stop("`", arg, "` must contain only numeric columns; non-numeric: ",
           paste(names(img)[!is_num], collapse = ", "), call. = FALSE)
    }
    img <- as.matrix(img)
  }
  if (is.vector(img) && is.numeric(img)) img <- matrix(img, ncol = 1)
  if (!is.matrix(img) || !is.numeric(img)) {
    stop("`", arg, "` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (nrow(img) == 0 || ncol(img) == 0) stop("`", arg, "` is empty.", call. = FALSE)
  if (anyNA(img) || any(!is.finite(img))) {
    stop("`", arg, "` contains missing or non-finite values. MIRTH currently ",
         "requires complete imaging data.", call. = FALSE)
  }
  if (is.null(colnames(img))) colnames(img) <- paste0("V", seq_len(ncol(img)))
  if (anyDuplicated(colnames(img))) {
    stop("`", arg, "` has duplicated column names.", call. = FALSE)
  }
  rownames(img) <- NULL
  img
}

.prepare_bat <- function(bat, n, arg) {
  if (is.data.frame(bat)) {
    if (ncol(bat) != 1) stop("`", arg, "` must be a vector.", call. = FALSE)
    bat <- bat[[1]]
  }
  if (length(bat) != n) {
    stop("`", arg, "` has length ", length(bat), " but the imaging data has ",
         n, " rows.", call. = FALSE)
  }
  if (anyNA(bat)) stop("`", arg, "` contains missing values.", call. = FALSE)
  droplevels(factor(bat))
}

.prepare_outcome <- function(y, n) {
  if (is.data.frame(y)) {
    if (ncol(y) != 1) stop("`train_outcome` must be a vector.", call. = FALSE)
    y <- y[[1]]
  }
  if (length(y) != n) {
    stop("`train_outcome` has length ", length(y), " but `train_img` has ", n,
         " rows.", call. = FALSE)
  }
  if (anyNA(y)) {
    stop("`train_outcome` contains missing values; drop these participants ",
         "from the training data.", call. = FALSE)
  }
  if (is.numeric(y)) {
    if (length(unique(y)) <= 2) {
      warning("`train_outcome` is numeric with <= 2 unique values and will be ",
              "treated as continuous. Convert it to a factor for classification.",
              call. = FALSE)
    }
    return(list(value = as.numeric(y), type = "regression", levels = NULL))
  }
  if (is.character(y) || is.logical(y)) y <- factor(y)
  if (!is.factor(y)) {
    stop("`train_outcome` must be numeric, factor, character, or logical.",
         call. = FALSE)
  }
  y <- droplevels(y)
  if (nlevels(y) < 2) stop("`train_outcome` has fewer than two classes.", call. = FALSE)
  if (is.ordered(y)) y <- factor(y, ordered = FALSE)
  list(value = y, type = "classification", levels = levels(y))
}

.prepare_covar <- function(covar, n, arg) {
  if (is.null(covar)) return(NULL)
  if (is.matrix(covar)) covar <- as.data.frame(covar)
  if (!is.data.frame(covar)) {
    stop("`", arg, "` must be a data frame (or NULL).", call. = FALSE)
  }
  covar <- as.data.frame(covar)  # drops tibble/data.table classes
  if (ncol(covar) == 0) return(NULL)
  if (nrow(covar) != n) {
    stop("`", arg, "` has ", nrow(covar), " rows but the imaging data has ",
         n, " rows.", call. = FALSE)
  }
  nm <- names(covar)
  bad <- nm[make.names(nm) != nm]
  if (length(bad)) {
    stop("Covariate names must be syntactic R names; rename: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(nm)) stop("`", arg, "` has duplicated column names.", call. = FALSE)
  clash <- intersect(nm, .reserved_names)
  if (length(clash)) {
    stop("Covariate names ", paste0("'", clash, "'", collapse = ", "),
         " are reserved by MIRTH/ComBat; please rename them.", call. = FALSE)
  }
  for (j in seq_along(covar)) {
    x <- covar[[j]]
    if (is.character(x) || is.logical(x)) covar[[j]] <- factor(x)
    else if (is.ordered(x)) covar[[j]] <- factor(x, ordered = FALSE)
    else if (!is.numeric(x) && !is.factor(x)) {
      stop("Covariate '", nm[j], "' must be numeric, factor, character, or logical.",
           call. = FALSE)
    }
  }
  rownames(covar) <- NULL
  covar
}

# Align new covariates to the training data by arranging columns in the same
# order and creating the same factor levels.
.align_covar <- function(new_covar, template, n) {
  if (is.null(template)) {
    if (!is.null(new_covar)) {
      stop("`new_covar` was supplied but the model was trained without covariates.",
           call. = FALSE)
    }
    return(NULL)
  }
  if (is.null(new_covar)) {
    message("`new_covar` not supplied: all covariates will be imputed.")
    new_covar <- template[rep(NA_integer_, n), , drop = FALSE]
    rownames(new_covar) <- NULL
    return(new_covar)
  }
  if (is.matrix(new_covar)) new_covar <- as.data.frame(new_covar)
  new_covar <- as.data.frame(new_covar)
  if (nrow(new_covar) != n) {
    stop("`new_covar` has ", nrow(new_covar), " rows but `new_img` has ", n,
         " rows.", call. = FALSE)
  }
  missing_cols <- setdiff(names(template), names(new_covar))
  if (length(missing_cols)) {
    stop("`new_covar` is missing column(s): ", paste(missing_cols, collapse = ", "),
         ". Supply them as NA to have them imputed.", call. = FALSE)
  }
  new_covar <- new_covar[, names(template), drop = FALSE]
  for (nm in names(template)) {
    x <- new_covar[[nm]]
    if (is.factor(template[[nm]])) {
      lev <- levels(template[[nm]])
      xc <- as.character(x)
      unknown <- setdiff(unique(xc[!is.na(xc)]), lev)
      if (length(unknown)) {
        stop("Covariate '", nm, "' has level(s) not seen in training: ",
             paste(unknown, collapse = ", "), call. = FALSE)
      }
      new_covar[[nm]] <- factor(xc, levels = lev)
    } else {
      if (!(is.numeric(x) || all(is.na(x)))) {
        stop("Covariate '", nm, "' must be numeric, as in the training data.",
             call. = FALSE)
      }
      new_covar[[nm]] <- as.numeric(x)
    }
  }
  rownames(new_covar) <- NULL
  new_covar
}

# Match new imaging columns to the training columns 
.align_img <- function(new_img, train_names) {
  new_img <- .prepare_img(new_img, "new_img")
  if (ncol(new_img) != length(train_names)) {
    stop("`new_img` has ", ncol(new_img), " columns but the model was trained on ",
         length(train_names), ".", call. = FALSE)
  }
  if (setequal(colnames(new_img), train_names)) {
    new_img <- new_img[, train_names, drop = FALSE]
  } else {
    warning("Column names of `new_img` do not match the training data; ",
            "matching columns by position.", call. = FALSE)
    colnames(new_img) <- train_names
  }
  new_img
}

.check_count <- function(x, arg, min = 1) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < min || x != round(x)) {
    stop("`", arg, "` must be a single whole number >= ", min, ".", call. = FALSE)
  }
  as.integer(x)
}
