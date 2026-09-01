## Shared resolution of covariate reference values for the setup functions.
##
## `contRef` and `catRef` both accept either a bare scalar, which applies to
## every covariate, or a named list with an entry per covariate plus a `default`
## component. The vocabulary is:
##
##   continuous   a number, "mean", "median", or "model"
##   categorical  a level, "mode", "lowest", or "model"
##
## "model" reads the reference out of the NONMEM control stream, using the same
## derivation as `createParamFunction()`, so a reference row built this way
## cannot disagree with the parameter function generated from the same model.

#' Most frequent value
#' @noRd
refMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Non-missing values of one covariate on deduplicated data
#' @noRd
refValues <- function(dedup, cov, missVal) {
  if (!cov %in% names(dedup)) {
    stop("Covariate ", cov, " is not present in the data.", call. = FALSE)
  }
  v <- dedup[[cov]]
  v <- v[v != missVal & !is.na(v)]
  if (length(v) == 0) stop("Covariate ", cov, " contains only missing values.",
                           call. = FALSE)
  v
}

#' Classify a covariate the same way `getCovStats()` does
#' @noRd
refCovType <- function(v, minLevels) {
  n <- length(unique(v))
  if (n > minLevels) "continuous" else if (n == 2) "binary" else
    if (n > 2) "multi" else "single"
}

#' Pull the spec for one covariate out of a scalar-or-list argument
#'
#' A bare scalar is itself the default. In a list, a named entry wins, then a
#' `default` component, then the caller's built-in fallback.
#'
#' @noRd
refSpec <- function(arg, cov, fallback) {
  if (is.null(arg)) return(fallback)
  if (!is.list(arg)) {
    if (length(arg) != 1) {
      stop("A `contRef` / `catRef` given as a vector must be a single value; ",
           "use a named list for per-covariate settings.", call. = FALSE)
    }
    return(arg)
  }
  if (!is.null(arg[[cov]])) return(arg[[cov]])
  if (!is.null(arg[["default"]])) return(arg[["default"]])
  fallback
}

#' Reference values derived from a NONMEM control stream
#'
#' `model` is either the list returned by [createParamFunction()] or a path to a
#' control stream, in which case the `$PK` block is parsed here with the same
#' rules. Returns a named list of `list(value, source, confident)`.
#'
#' @noRd
refModelValues <- function(model, missVal) {
  if (is.null(model)) {
    stop("`model` is required when a reference is set to \"model\". Supply the ",
         "control stream path, or the list returned by createParamFunction().",
         call. = FALSE)
  }
  if (is.list(model)) {
    if (is.null(model$covRef)) {
      stop("`model` must be a control stream path or the list returned by ",
           "createParamFunction().", call. = FALSE)
    }
    return(model$covRef)
  }
  if (!is.character(model) || length(model) != 1) {
    stop("`model` must be a control stream path or the list returned by ",
         "createParamFunction().", call. = FALSE)
  }

  mod <- nmReadModel(model)
  pk  <- nmRecord(mod, "\\$PK\\b")
  if (nrow(pk) == 0) {
    stop("No $PK record found in ", basename(model),
         ", so no reference values can be read from it.", call. = FALSE)
  }
  stmts <- nmSimplifyStmts(nmParseStatements(pk, model))
  syms  <- nmSymbols(stmts)
  covs  <- intersect(syms$used, setdiff(nmInputNames(mod), syms$assigned))
  nmCovRef(stmts, covs, missVal)
}

#' Resolve a reference value for every covariate
#'
#' Returns a named list of `list(value, source, confident)`, one entry per
#' covariate. `catFallback` is the built-in default for categorical covariates:
#' `"mode"` where a reference subject is wanted, `"lowest"` where the question is
#' which level to drop when one-hot encoding.
#'
#' @noRd
refResolve <- function(data, covariates, contRef = "median", catRef = NULL,
                       model = NULL, minLevels = 10, idVar = "ID",
                       missVal = -99, nsig = 3, catFallback = "mode") {

  dedup <- data %>% dplyr::distinct(!!rlang::sym(idVar), .keep_all = TRUE)

  ## Parse the control stream once, and only if something actually asks for it.
  wants <- function(arg) {
    if (is.null(arg)) return(FALSE)
    any(vapply(if (is.list(arg)) arg else list(arg),
               function(x) is.character(x) && identical(x[1], "model"),
               logical(1)))
  }
  modelRefs <- if (wants(contRef) || wants(catRef)) {
    refModelValues(model, missVal)
  } else NULL

  out <- list()
  for (cov in covariates) {
    v    <- refValues(dedup, cov, missVal)
    type <- refCovType(v, minLevels)
    cont <- type == "continuous"
    spec <- refSpec(if (cont) contRef else catRef, cov,
                    if (cont) "median" else catFallback)

    if (is.numeric(spec)) {
      # An explicit value is used as given, never rounded.
      out[[cov]] <- list(value = spec, source = "supplied directly",
                         confident = TRUE)
      next
    }
    if (!is.character(spec) || length(spec) != 1) {
      stop("The reference setting for covariate ", cov,
           " must be a single value or one of \"mean\", \"median\", \"mode\", ",
           "\"lowest\" or \"model\".", call. = FALSE)
    }

    out[[cov]] <- switch(spec,
      "model" = {
        r <- modelRefs[[cov]]
        if (is.null(r)) {
          stop("No reference value for covariate ", cov,
               " could be derived from the model. Give it explicitly, e.g. ",
               if (cont) paste0("contRef = list(", cov, " = <value>)")
               else paste0("catRef = list(", cov, " = <level>)"), ".",
               call. = FALSE)
        }
        r
      },
      "mean" = {
        if (!cont) stop("\"mean\" is not a reference for the categorical ",
                        "covariate ", cov, ".", call. = FALSE)
        list(value = signif(mean(v), nsig), source = "mean of the data",
             confident = TRUE)
      },
      "median" = {
        if (!cont) stop("\"median\" is not a reference for the categorical ",
                        "covariate ", cov, ".", call. = FALSE)
        list(value = signif(stats::median(v), nsig),
             source = "median of the data", confident = TRUE)
      },
      "mode" = {
        if (cont) stop("\"mode\" is not a reference for the continuous ",
                       "covariate ", cov, ".", call. = FALSE)
        list(value = refMode(v), source = "most common level in the data",
             confident = TRUE)
      },
      "lowest" = {
        if (cont) stop("\"lowest\" is not a reference for the continuous ",
                       "covariate ", cov, ".", call. = FALSE)
        list(value = sort(unique(v))[1], source = "lowest level in the data",
             confident = TRUE)
      },
      stop("Unknown reference setting \"", spec, "\" for covariate ", cov,
           ". Use a value, \"mean\", \"median\", \"mode\", \"lowest\" or ",
           "\"model\".", call. = FALSE)
    )
  }
  out
}

#' Warn about the deprecated `refLevels` argument and fold it into `catRef`
#' @noRd
refLevelsToCatRef <- function(refLevels, catRef, fn) {
  if (is.null(refLevels)) return(catRef)
  if (!is.null(catRef)) {
    stop("Supply either `catRef` or the deprecated `refLevels`, not both.",
         call. = FALSE)
  }
  warning("`refLevels` is deprecated in ", fn, "(); use `catRef` instead. ",
          "`catRef` additionally accepts \"mode\", \"lowest\", \"model\" and a ",
          "`default` component.", call. = FALSE)
  refLevels
}

#' Reference level used when one-hot encoding a multi-level categorical
#'
#' Resolves `catRef` for the encoding question specifically - which level is
#' dropped - so the fallback is the lowest level rather than the mode. Returns
#' the level itself.
#'
#' @noRd
refEncodingLevel <- function(catRef, cov, levs, model, missVal) {
  spec <- refSpec(catRef, cov, "lowest")
  lev <- if (is.numeric(spec)) {
    spec
  } else if (identical(spec, "lowest")) {
    levs[1]
  } else if (identical(spec, "mode")) {
    NULL                       # resolved by the caller, which has the data
  } else if (identical(spec, "model")) {
    r <- refModelValues(model, missVal)[[cov]]
    if (is.null(r)) {
      stop("No reference level for covariate ", cov,
           " could be derived from the model. Give it explicitly, e.g. ",
           "catRef = list(", cov, " = <level>).", call. = FALSE)
    }
    r$value
  } else {
    stop("Unknown reference setting \"", spec, "\" for covariate ", cov, ".",
         call. = FALSE)
  }
  lev
}
