#' Setup Covariate Data Frame for Forest Plots
#'
#' @description A high-level wrapper that streamlines the creation of the covariate
#'   input data frame for `getForestDF` functions. It calculates summary statistics
#'   and reshapes the output into the required plotting format. It natively supports
#'   `conditionalCovs` (e.g., for FREM workflows) and allows explicit control over
#'   how inactive covariate cells are populated.
#'
#' @details
#'   **Methodological Note on Reference Calculations:** To maintain mathematical
#'   symmetry with how primary `covariates` are evaluated, the reference values
#'   for `conditionalCovs` (and primary covariates when `useMissVal = FALSE`) are
#'   calculated strictly on **deduplicated baseline data** (one record per `idVar`).
#'   This prevents subjects with dense longitudinal sampling from skewing the reference
#'   values. If you require longitudinal aggregation, you must pre-process your
#'   dataset before passing it to this function or post-process the output data.frame.
#'
#' @param data A data frame that includes the covariates to summarize. Only the
#'   first record per subject (identified by `idVar`) will be used for both
#'   primary quantile derivation and reference calculations.
#' @param covariates A character vector of primary covariate names to evaluate.
#' @param conditionalCovs A character vector of supplementary covariates. These are
#'   expanded into their own rows like primary covariates, but when "inactive",
#'   they default to their dataset reference value (mode for categorical, mean/median
#'   for continuous) rather than `missVal`. Defaults to `NULL`.
#' @param useMissVal Logical. If `TRUE` (default), inactive primary `covariates`
#'   retain `missVal`. If `FALSE`, inactive primary covariates are replaced with
#'   their computed baseline references.
#' @param contRef Reference setting for continuous covariates: a number,
#'   `"mean"`, `"median"` (default) or `"model"`, or a named list of those with
#'   an optional `default` component, e.g.
#'   `list(WT = 75, AGE = "mean", default = "median")`. Applies to
#'   `conditionalCovs`, and to the primary covariates when `useMissVal = FALSE`.
#' @param minLevels The maximum number of unique values a covariate can have to be
#'   treated as categorical. Default is 10.
#' @param probs A numeric vector of two probabilities used to calculate quantiles
#'   for continuous covariates. Defaults to `c(0.05, 0.95)`.
#' @param idVar The name of the subject identifier column. Defaults to `"ID"`.
#' @param missVal The numeric value indicating missing data, used to fill inactive
#'   states. Defaults to -99.
#' @param nsig The number of significant digits for rounding continuous covariates.
#'   Defaults to 3.
#' @param catRef Reference setting for categorical covariates: a level,
#'   `"mode"`, `"lowest"` or `"model"`, or a named list of those with an optional
#'   `default` component, e.g. `list(GENO = 2)`. It selects both the level
#'   represented by the all-zero one-hot row and, for background cells, the
#'   reference level. Pass the same `catRef` to [setupDfRefRow()] so the columns
#'   line up. Replaces `refLevels`.
#' @param model The NONMEM control stream, required when a reference is set to
#'   `"model"`. Either a path to the `.mod` file or the list returned by
#'   [createParamFunction()].
#' @param refLevels Deprecated. Use `catRef`.
#' @param sep The separator between the covariate name and the level in the
#'   one-hot column names. Defaults to `"_"`.
#'
#' @return A data frame formatted for use in empirical and SCM forest plot functions.
#' @export
#'
#' @seealso \code{\link{getCovStats}}, \code{\link{createInputForestData}}
#'   `vignette("Part3-deep-dive-forest-plot-inputs", package = "PMXForest")` for how this fits the whole
#'   workflow.
#'
#' @examples
#' dfData <- read.csv(
#'   system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
#' )
#'
#' # Continuous covariates -> 5th/95th percentiles; GENO -> one row per level.
#' # Inactive cells hold missVal (-99); the parameter function handles that.
#' setupDfCovs(dfData, covariates = c("WT", "AGE", "SEX", "GENO"), idVar = "ID")
#'
#' # CRCL as an conditionalCov: it sits at its median on the WT/SEX rows
#' # instead of at missVal.
#' setupDfCovs(dfData,
#'   covariates = c("WT", "SEX"), conditionalCovs = "CRCL",
#'   idVar = "ID"
#' )
#'
#' # useMissVal = FALSE: inactive primary covariates also hold their reference
#' setupDfCovs(dfData,
#'   covariates = c("WT", "SEX"), useMissVal = FALSE,
#'   idVar = "ID"
#' )
#'
#' # Align the GENO one-hot columns with a model whose reference genotype is 2
#' setupDfCovs(dfData,
#'   covariates = c("WT", "GENO"), catRef = list(GENO = 2),
#'   idVar = "ID"
#' )
setupDfCovs <- function(data, covariates, conditionalCovs = NULL, useMissVal = TRUE,
                        contRef = "median", catRef = NULL, model = NULL,
                        refLevels = NULL, minLevels = 10,
                        probs = c(0.05, 0.95), idVar = "ID",
                        missVal = -99, nsig = 3, sep = "_") {
  catRef <- refLevelsToCatRef(refLevels, catRef, "setupDfCovs")
  all_covs <- unique(c(covariates, conditionalCovs))

  # 1. Calculate statistical summaries
  stats_list <- getCovStats(
    data = data,
    covariates = all_covs,
    minLevels = minLevels,
    probs = probs,
    idVar = idVar,
    missVal = missVal,
    nsig = nsig,
    catRef = catRef,
    model = model,
    sep = sep
  )

  # 2. Reshape into the basic forest plot data.frame structure
  df_covs <- createInputForestData(
    listCovs = stats_list,
    iMiss = missVal
  )

  # 3. Determine which covariates require background reference replacement
  target_covs <- NULL
  if (!is.null(conditionalCovs)) {
    target_covs <- setdiff(conditionalCovs, covariates)
  }
  if (!useMissVal) {
    target_covs <- unique(c(target_covs, covariates))
  }

  # 4. Post-process to replace background missVal with computed references
  if (length(target_covs) > 0) {
    refs <- refResolve(data, target_covs,
      contRef = contRef, catRef = catRef,
      model = model, minLevels = minLevels, idVar = idVar,
      missVal = missVal, nsig = nsig, catFallback = "mode"
    )

    dedup_data <- data %>% dplyr::distinct(!!rlang::sym(idVar), .keep_all = TRUE)

    for (acov in target_covs) {
      v <- refValues(dedup_data, acov, missVal)
      type <- refCovType(v, minLevels)
      value <- refs[[acov]]$value

      if (type == "multi") {
        levs <- sort(unique(v))
        encLev <- refEncodingLevel(catRef, acov, levs, model, missVal)
        if (is.null(encLev)) encLev <- refMode(v)
        for (lev in setdiff(levs, encLev)) {
          col_name <- paste0(acov, sep, lev)
          if (col_name %in% names(df_covs)) {
            df_covs[[col_name]][df_covs[[col_name]] == missVal] <-
              as.numeric(value == lev)
          }
        }
      } else {
        df_covs[[acov]][df_covs[[acov]] == missVal] <- value
      }
    }
  }

  return(df_covs)
}
