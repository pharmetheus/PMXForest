#' One-Hot Encode Categorical Covariates
#'
#' @description Adds one-hot (dummy) columns for multi-level categorical
#'   covariates to a data frame, following the naming convention used across the
#'   Pharmetheus forest-plot tooling: `<covariate><sep><level>`, with the lowest
#'   level dropped as the reference unless stated otherwise.
#'
#'   It is the shared building block behind `setupDfCovs()` and the optional
#'   `oneHot` argument of `getForestDFSCM()` and `getForestDFemp()`, and it
#'   matches the columns produced by NONMEM FREM data sets and
#'   `PMXFrem::createFREMModel()`.
#'
#' @details
#'   For a covariate `GENO` with non-missing levels `1, 2, 3, 4` and the default
#'   reference (level 1), the function adds `GENO_2`, `GENO_3`, and `GENO_4`; the
#'   reference level is represented by all three being `0`.
#'
#'   The function is idempotent: if a target column already exists it is left
#'   untouched when its non-missing rows agree with the requested encoding, and
#'   an error is raised when they do not. This makes it safe to call on data sets
#'   that already carry FREM dummy columns.
#'
#' @param data A data frame containing the raw categorical columns.
#' @param spec The covariates to encode, in one of two forms:
#'   \itemize{
#'     \item A character vector of column names. Each is encoded with its lowest
#'       non-missing level as the reference.
#'     \item A named list. Each element is either `NULL` (lowest level as
#'       reference), a single value (used as the reference level), or a list with
#'       `ref` (reference level) and optional `levels` (the non-reference levels
#'       to create columns for; defaults to all remaining levels).
#'   }
#' @param sep The separator between the covariate name and the level in the new
#'   column names. Defaults to `"_"`, matching the FREM convention. Changing it
#'   breaks alignment with FREM tooling and is discouraged.
#' @param missVal The value indicating missing data. Rows holding `missVal` (or
#'   `NA`) in the raw column receive `missVal` in the dummy columns unless
#'   `imputeMissing = TRUE`. Defaults to -99.
#' @param includeReference Logical. If `TRUE`, a dummy column is also created for
#'   the reference level. Defaults to `FALSE`.
#' @param imputeMissing Logical. If `TRUE`, rows with a missing raw value are set
#'   to `0` in every dummy column (imputed to the reference category). If `FALSE`
#'   (default), they are set to `missVal`, so downstream code can keep treating
#'   them as missing.
#' @param dropOriginal Logical. If `TRUE`, the raw categorical column is removed
#'   after its dummy columns have been added. Defaults to `FALSE`.
#'
#' @return `data` with the dummy columns added (numeric `0`/`1`, or `missVal` on
#'   missing rows when `imputeMissing = FALSE`). If a `COVARIATEGROUPS` column is
#'   present it is moved back to the last position.
#' @export
#'
#' @seealso \code{\link{setupDfCovs}}, \code{\link{getCovStats}}
#'
#' @examples
#' df <- data.frame(
#'   ID   = 1:6,
#'   GENO = c(1, 2, 3, 4, 2, -99),
#'   RACE = c(1, 1, 2, 3, 2, 1)
#' )
#'
#' # Default: lowest level is the reference, missing rows kept as missVal
#' oneHotEncode(df, spec = c("GENO", "RACE"))
#'
#' # Explicit reference level for GENO (level 2), and drop the raw column
#' oneHotEncode(df, spec = list(GENO = list(ref = 2)), dropOriginal = TRUE)
#'
#' # Impute missing GENO to the reference category
#' oneHotEncode(df, spec = "GENO", imputeMissing = TRUE)
oneHotEncode <- function(data, spec, sep = "_", missVal = -99,
                         includeReference = FALSE, imputeMissing = FALSE,
                         dropOriginal = FALSE) {

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")

  normSpec <- normalizeOneHotSpec(spec, data, missVal, includeReference)

  for (cov in names(normSpec)) {
    ref     <- normSpec[[cov]]$ref
    raw     <- data[[cov]]
    isMiss  <- is.na(raw) | raw == missVal

    for (lev in normSpec[[cov]]$levels) {
      col <- paste0(cov, sep, lev)

      newVal <- if (imputeMissing) {
        as.numeric(!isMiss & raw == lev)
      } else {
        ifelse(isMiss, missVal, as.numeric(raw == lev))
      }

      if (col %in% names(data)) {
        ## Idempotency: compare only on rows that are not missing in the raw
        ## column, so a pre-existing imputed (missing -> 0) encoding is accepted.
        keep <- if (imputeMissing) seq_along(newVal) else which(!isMiss)
        ok <- isTRUE(all.equal(as.numeric(data[[col]][keep]),
                               as.numeric(newVal[keep])))
        if (!ok) {
          stop("Column '", col, "' already exists and is inconsistent with the ",
               "requested one-hot encoding of '", cov, "'.")
        }
      } else {
        data[[col]] <- newVal
      }
    }

    if (dropOriginal) data[[cov]] <- NULL
  }

  if ("COVARIATEGROUPS" %in% names(data)) {
    data <- data[c(setdiff(names(data), "COVARIATEGROUPS"), "COVARIATEGROUPS")]
  }

  data
}

#' Normalise a one-hot encoding specification
#'
#' Internal helper. Turns the user-facing `spec` argument of [oneHotEncode()]
#' into a named list of `list(ref = , levels = )` entries, dropping covariates
#' that are absent or have fewer than two non-missing levels (with a warning).
#'
#' @inheritParams oneHotEncode
#' @return A named list, one entry per usable covariate.
#' @keywords internal
#' @noRd
normalizeOneHotSpec <- function(spec, data, missVal, includeReference) {

  if (is.character(spec)) {
    covs <- spec
    spec <- vector("list", length(covs))
    names(spec) <- covs
  }

  if (!is.list(spec) || is.null(names(spec)) || any(names(spec) == "")) {
    stop("`spec` must be a character vector of covariate names or a named list.")
  }

  out <- list()
  for (cov in names(spec)) {

    if (!cov %in% names(data)) {
      warning("One-hot encoding: covariate '", cov, "' not found in the data; skipped.")
      next
    }

    obs  <- data[[cov]]
    obs  <- obs[!is.na(obs) & obs != missVal]
    levs <- sort(unique(obs))

    if (length(levs) < 2) {
      warning("One-hot encoding: covariate '", cov,
              "' has fewer than two non-missing levels; skipped.")
      next
    }

    el <- spec[[cov]]
    wantLevels <- NULL

    if (is.null(el)) {
      ref <- levs[1]
    } else if (is.list(el)) {
      ref <- if (!is.null(el$ref)) el$ref else levs[1]
      wantLevels <- el$levels
    } else if (length(el) == 1) {
      ref <- el
    } else {
      stop("One-hot spec for '", cov,
           "' must be NULL, a single reference value, or a list with `ref`/`levels`.")
    }

    if (!ref %in% levs) {
      warning("One-hot encoding: reference level ", ref, " for covariate '", cov,
              "' is not present in the data.")
    }

    if (is.null(wantLevels)) {
      wantLevels <- if (includeReference) levs else setdiff(levs, ref)
    } else if (includeReference && !ref %in% wantLevels) {
      wantLevels <- c(ref, wantLevels)
    }

    out[[cov]] <- list(ref = ref, levels = wantLevels)
  }

  out
}
