#' Get Covariate Statistics for Forest Plots
#'
#' @description Quickly generates a list of summary statistics for specified
#'   covariates, formatted for use in forest plots. The function handles
#'   continuous, binary, and multi-level categorical variables differently based
#'   on their number of unique values.
#'
#' @param data A data frame that includes the covariates to summarise. Only the
#'   first record per subject (identified by `idVar`) will be used.
#' @param covariates A character vector of covariate names as they appear in the
#'   data frame.
#' @param minLevels The maximum number of unique values a covariate can have to be
#'   treated as categorical. Covariates with more unique values than `minLevels`
#'   are treated as continuous. Default is 10.
#' @param probs A numeric vector of two probabilities used to calculate quantiles
#'   for continuous covariates. Defaults to `c(0.05, 0.95)`.
#' @param idVar The name of the subject identifier column, used to select one
#'   record per subject. Defaults to `"ID"`.
#' @param missVal The value indicating missing data, which will be excluded from
#'   calculations. Defaults to -99.
#' @param nsig The number of significant digits (passed to the `signif`
#'   function) for rounding the summary of continuous covariates. Defaults to 3.
#' @param refLevels An optional named list (or named vector) giving the reference
#'   level for individual multi-level categorical covariates, e.g.
#'   `list(GENO = 2)`. Covariates not named here use their lowest level as the
#'   reference. Has no effect on continuous or binary covariates.
#' @param sep The separator between the covariate name and the level in the
#'   one-hot column names for multi-level categorical covariates. Defaults to
#'   `"_"`, matching the NONMEM FREM convention.
#'
#' @return A list where each element corresponds to a covariate.
#'   \itemize{
#'     \item For **continuous** covariates, the element is a named numeric vector
#'       of quantiles.
#'     \item For **binary** covariates, the element is a sorted vector of the two
#'       unique values.
#'     \item For **multi-level categorical** covariates, the element is a nested
#'       list of one-hot encoded vectors. The reference level (the lowest level,
#'       or the one named in `refLevels`) is represented by all vectors being 0.
#'   }
#' @export
#'
#' @examples
#' # --- Basic Example ---
#'
#' # 1. Create a sample dataset
#' # This data includes duplicate IDs, continuous (WT), binary (SEX),
#' # multi-level categorical (RACE), and a variable with a missing value (BMI).
#' sample_data <- data.frame(
#'   ID = c(1, 1, 2, 3, 4, 5),
#'   WT = c(60.5, 61.0, 70.2, 80.8, 65.1, 90.3),
#'   SEX = c(0, 0, 1, 0, 1, 0),
#'   RACE = c(1, 1, 2, 3, 1, 2),
#'   BMI = c(22.1, 22.4, 25.3, -99, 23.5, 28.9)
#' )
#'
#' # 2. Define covariates and get statistics
#' covariates_to_summarize <- c("WT", "SEX", "RACE")
#'
#' # Note: WT has 5 unique values in the distinct-ID dataset.
#' # To ensure it's treated as continuous, we set minLevels to be less than 5.
#' cov_stats <- getCovStats(
#'   data = sample_data,
#'   covariates = covariates_to_summarize,
#'   minLevels = 4
#' )
#'
#' # 3. View the output list structure
#' print(cov_stats)
getCovStats <- function (data, covariates, minLevels = 10, probs = c(0.05, 0.95),
                         idVar = "ID", missVal = -99, nsig = 3,
                         refLevels = NULL, sep = "_") {

  # REFACTORED: Use sym() instead of ensym() to allow programmatic wrapping.
  data <- data %>% distinct(!!sym(idVar), .keep_all = TRUE)

  ## Check the input
  if(!all(covariates %in% names(data))) stop("Not all covariates are present in the data.")

  retList <- list()
  for (myCov in covariates) {
    dataTmp <- data[data[[myCov]] != missVal, ]
    if (length(unique(dataTmp[[myCov]])) <= minLevels) {
      numLevs <- length(unique(dataTmp[[myCov]]))
      if (numLevs == 2) {
        # REFACTORED: Sort binary variables for predictable output ordering
        retList[[myCov]] <- sort(unique(dataTmp[[myCov]]))
      }
      else {
        levs    <- sort(unique(dataTmp[[myCov]]))
        numLevs <- length(levs)
        if (numLevs < 2) {
          ## Degenerate: fewer than two non-missing levels, no contrast to form.
          retList[[myCov]] <- levs
        } else {
          refLev <- if (!is.null(refLevels[[myCov]])) refLevels[[myCov]] else levs[1]
          if (!refLev %in% levs) {
            stop("Reference level ", refLev, " for covariate '", myCov,
                 "' is not present in the data.")
          }
          covList <- list()
          for (i in seq_len(numLevs)) {
            if (levs[i] == refLev) next
            vec <- rep(0, numLevs)
            vec[i] <- 1
            covList[[paste0(myCov, sep, levs[i])]] <- vec
          }
          retList[[myCov]] <- covList
        }
      }
    }
    else {
      retList[[myCov]] <- signif(quantile(dataTmp[[myCov]], p=probs), digits = nsig)
    }
  }
  return(retList)
}
