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
#' @param catRef The reference level for multi-level categorical covariates -
#'   the level represented by the all-zero one-hot row. Either a single setting
#'   applied to every covariate or a named list with an optional `default`
#'   component; each entry is a level, `"lowest"` (the default), `"mode"` or
#'   `"model"`. For example `list(GENO = 2)`. Has no effect on continuous or
#'   binary covariates.
#' @param model The NONMEM control stream, required when `catRef` is `"model"`.
#'   Either a path to the `.mod` file or the list returned by
#'   [createParamFunction()].
#' @param refLevels Deprecated. Use `catRef`.
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
#'       or the one named in `catRef`) is represented by all vectors being 0.
#'   }
#' @export
#'
#' @examples
#' dfData <- read.csv(
#'   system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
#' )
#'
#' # Continuous covariates -> quantiles; binary -> the two levels;
#' # GENO (4 levels) -> a one-hot list with level 1 as the reference.
#' getCovStats(dfData, covariates = c("WT", "AGE", "SEX", "GENO"), idVar = "ID")
#'
#' # Non-default quantiles and rounding for a continuous covariate
#' getCovStats(dfData, "CRCL", idVar = "ID", probs = c(0.1, 0.9), nsig = 4)
#'
#' # Use genotype level 2 as the reference instead of the lowest level
#' getCovStats(dfData, "GENO", idVar = "ID", catRef = list(GENO = 2))
getCovStats <- function (data, covariates, minLevels = 10, probs = c(0.05, 0.95),
                         idVar = "ID", missVal = -99, nsig = 3,
                         catRef = NULL, model = NULL, refLevels = NULL,
                         sep = "_") {

  catRef <- refLevelsToCatRef(refLevels, catRef, "getCovStats")

  # REFACTORED: Use sym() instead of ensym() to allow programmatic wrapping.
  data <- data %>% distinct(!!sym(idVar), .keep_all = TRUE)

  ## Check the input
  if(!all(covariates %in% names(data))) stop("Not all covariates are present in the data.")

  retList <- list()
  for (myCov in covariates) {
    ## `x != missVal` is NA where x is NA, and indexing rows by a logical NA
    ## keeps an all-NA row rather than dropping it. That leaked NA used to be
    ## counted as a level (turning a binary covariate into a one-hot list) and
    ## to reach quantile(), whose na.rm is FALSE. Drop it explicitly, as
    ## refValues() and setupCovExpressionsList() already do.
    dataTmp <- data[!is.na(data[[myCov]]) & data[[myCov]] != missVal, ]
    if (nrow(dataTmp) == 0) {
      stop("Covariate ", myCov, " contains only missing values.", call. = FALSE)
    }
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
          refLev <- refEncodingLevel(catRef, myCov, levs, model, missVal)
          if (is.null(refLev)) refLev <- refMode(dataTmp[[myCov]])
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
