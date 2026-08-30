#' Setup Covariate Data Frame for Forest Plots
#'
#' @description A high-level wrapper that streamlines the creation of the covariate
#'   input data frame for `getForestDF` functions. It calculates summary statistics
#'   and reshapes the output into the required plotting format. It natively supports
#'   `additionalCovs` (e.g., for FREM workflows) and allows explicit control over
#'   how inactive covariate cells are populated.
#'
#' @details
#'   **Methodological Note on Reference Calculations:** To maintain mathematical
#'   symmetry with how primary `covariates` are evaluated, the reference values
#'   for `additionalCovs` (and primary covariates when `useMissVal = FALSE`) are
#'   calculated strictly on **deduplicated baseline data** (one record per `idVar`).
#'   This prevents subjects with dense longitudinal sampling from skewing the reference
#'   values. If you require longitudinal aggregation, you must pre-process your
#'   dataset before passing it to this function or post-process the output data.frame.
#'
#' @param data A data frame that includes the covariates to summarize. Only the
#'   first record per subject (identified by `idVar`) will be used for both
#'   primary quantile derivation and reference calculations.
#' @param covariates A character vector of primary covariate names to evaluate.
#' @param additionalCovs A character vector of supplementary covariates. These are
#'   expanded into their own rows like primary covariates, but when "inactive",
#'   they default to their dataset reference value (mode for categorical, mean/median
#'   for continuous) rather than `missVal`. Defaults to `NULL`.
#' @param useMissVal Logical. If `TRUE` (default), inactive primary `covariates`
#'   retain `missVal`. If `FALSE`, inactive primary covariates are replaced with
#'   their computed baseline references.
#' @param contRef A character string indicating how to calculate the reference
#'   value for continuous covariates. Must be either `"median"` or `"mean"`.
#'   Defaults to `"median"`.
#' @param minLevels The maximum number of unique values a covariate can have to be
#'   treated as categorical. Default is 10.
#' @param probs A numeric vector of two probabilities used to calculate quantiles
#'   for continuous covariates. Defaults to `c(0.05, 0.95)`.
#' @param idVar The name of the subject identifier column. Defaults to `"ID"`.
#' @param missVal The numeric value indicating missing data, used to fill inactive
#'   states. Defaults to -99.
#' @param nsig The number of significant digits for rounding continuous covariates.
#'   Defaults to 3.
#' @param refLevels An optional named list giving the reference level for
#'   individual multi-level categorical covariates, e.g. `list(GENO = 2)`. Use
#'   this to align the one-hot columns with a model whose reference genotype (or
#'   race, etc.) is not the lowest level. Covariates not named here use their
#'   lowest level as the reference. Passed through to `getCovStats()`.
#' @param sep The separator between the covariate name and the level in the
#'   one-hot column names. Defaults to `"_"`.
#'
#' @return A data frame formatted for use in empirical and SCM forest plot functions.
#' @export
#'
#' @seealso \code{\link{getCovStats}}, \code{\link{createInputForestData}}
#'
#' @examples
#' # 1. Create a sample dataset
#' sample_data <- data.frame(
#'   ID = c(1, 2, 3, 4, 5, 6),
#'   WT = c(60.5, 70.2, 70.5, 80.8, 65.1, 90.0), # 6 unique values
#'   SEX = c(0, 1, 1, 0, 1, 0),                  # 2 levels
#'   RACE = c(1, 2, 2, 3, 1, 3),                 # 3 levels
#'   FOOD = c(1, 0, 0, 1, 1, 0)                  # 2 levels (Mode is 1)
#' )
#'
#' # 2. Generate the formatted covariate data frame
#' # WT, SEX, and RACE are primary covariates.
#' # FOOD is a supplementary additional covariate.
#' # Setting minLevels = 4 ensures WT (6 unique values) remains continuous,
#' # while RACE (3 unique values) is correctly routed to the multi-level pathway.
#'
#' # Default. Relying on the parameter function to handle missing covariates.
#' # Inactive primary covariates hold missVal (-99).
#' # The additional covariate (FOOD) holds its baseline mode (1) in the background.
#' df_def <- setupDfCovs(
#'   data = sample_data,
#'   covariates = c("WT", "SEX", "RACE"),
#'   additionalCovs = "FOOD",
#'   useMissVal = TRUE,
#'   minLevels = 4
#' )
#' print(df_def)
#'
#' # Alternative: All covariate values are handled by dfCovs.
#' # Inactive primary covariates are ALSO replaced with their computed baseline references.
#' # When WT is the active varying covariate, SEX sits at 0, RACE_2 at 0, RACE_3 at 0, etc.
#' df_alt <- setupDfCovs(
#'   data = sample_data,
#'   covariates = c("WT", "SEX", "RACE"),
#'   additionalCovs = "FOOD",
#'   useMissVal = FALSE,
#'   minLevels = 4
#' )
#' print(df_alt)
setupDfCovs <- function(data, covariates, additionalCovs = NULL, useMissVal = TRUE,
                        contRef = c("median", "mean"), minLevels = 10,
                        probs = c(0.05, 0.95), idVar = "ID",
                        missVal = -99, nsig = 3, refLevels = NULL, sep = "_") {

  contRef <- match.arg(contRef)
  all_covs <- unique(c(covariates, additionalCovs))

  # 1. Calculate statistical summaries
  stats_list <- getCovStats(
    data = data,
    covariates = all_covs,
    minLevels = minLevels,
    probs = probs,
    idVar = idVar,
    missVal = missVal,
    nsig = nsig,
    refLevels = refLevels,
    sep = sep
  )

  # 2. Reshape into the basic forest plot data.frame structure
  df_covs <- createInputForestData(
    listCovs = stats_list,
    iMiss = missVal
  )

  # 3. Determine which covariates require background reference replacement
  target_covs <- NULL
  if (!is.null(additionalCovs)) {
    target_covs <- setdiff(additionalCovs, covariates)
  }
  if (!useMissVal) {
    target_covs <- unique(c(target_covs, covariates))
  }

  # 4. Post-process to replace background missVal with computed references
  if (length(target_covs) > 0) {

    dedup_data <- data %>% dplyr::distinct(!!rlang::sym(idVar), .keep_all = TRUE)

    get_mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }

    for (acov in target_covs) {
      v <- dedup_data[[acov]][dedup_data[[acov]] != missVal & !is.na(dedup_data[[acov]])]

      if (length(v) == 0) stop(paste("Covariate", acov, "contains only missing values."))

      n_levs <- length(unique(v))

      if (n_levs <= minLevels) {
        mode_val <- get_mode(v)

        if (n_levs == 2) {
          df_covs[[acov]][df_covs[[acov]] == missVal] <- mode_val
        } else {
          levs   <- sort(unique(v))
          refLev <- if (!is.null(refLevels[[acov]])) refLevels[[acov]] else levs[1]
          for (lev in setdiff(levs, refLev)) {
            col_name <- paste0(acov, sep, lev)
            ref_val  <- ifelse(mode_val == lev, 1, 0)
            if (col_name %in% names(df_covs)) {
              df_covs[[col_name]][df_covs[[col_name]] == missVal] <- ref_val
            }
          }
        }
      } else {
        ref_val <- if (contRef == "median") median(v) else mean(v)
        ref_val <- signif(ref_val, nsig)
        df_covs[[acov]][df_covs[[acov]] == missVal] <- ref_val
      }
    }
  }

  return(df_covs)
}
