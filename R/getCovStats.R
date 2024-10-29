#' getCovStats
#'
#' @description To quickly get a dfCovs data from for Forest plots by just supplying a data file and a vector of covariate names.
#' @param data A data frame that includes the covariates to summarise. Only the first line per subject will be used in the summary.
#' @param covariates A character vector of covariate names as they appear in the data file.
#' @param minLevels The maximum number of unique values for a covariate to be regarded as continus. Devfault is 10.
#' @param probs A vector of two number - the lower and upper percentiles to use for summarising continuous covariates.
#' @param idVar The name of the ID column to use for identifying one record per subject.
#' @param missVal The missing value indicator. Default is -99.
#' @param nsig The number of `digits` passed to the `signif` function when applied
#' to covariates with more than `minLevels` of unique values.
#'
#' @return A data frame with covariate statistics to be used as an argument for Forest plots.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' data <- read.csv("inst/extdata/DAT-1-MI-PMX-2.csv")
#'
#' dfCovs <- createInputForestData(
#'   getCovStats(data,c("WT","BMI","SEX"),missVal=-99),
#'   iMiss    =-99
#' )}
#'
getCovStats <- function (data, covariates, minLevels = 10, probs = c(0.05, 0.95),
                         idVar = "ID", missVal = -99, nsig = 3) {

  data <- data %>% distinct(!!ensym(idVar), .keep_all = TRUE)

  ## Check the input
  if(!all(covariates %in% names(data))) stop("Not all covariates are present in the data.")

  data <- subset(data, !duplicated(idVar))
  retList <- list()
  for (myCov in covariates) {
    dataTmp <- data[data[[myCov]] != missVal, ]
    if (length(unique(dataTmp[[myCov]])) <= minLevels) {
      numLevs <- length(unique(dataTmp[[myCov]]))
      if (numLevs == 2) {
        retList[[myCov]] <- unique(dataTmp[[myCov]])
      }
      else {
        levs <- sort(unique(dataTmp[[myCov]]))
        covList <- list()
        for (i in 2:numLevs) {
          vec <- rep(0, numLevs)
          vec[i] <- 1
          covList[[paste0(myCov, "_", levs[i])]] <- vec
        }
        retList[[myCov]] <- covList
      }
    }
    else {
      retList[[myCov]] <- signif(quantile(dataTmp[[myCov]],p=probs),digits = nsig)
    }
  }
  return(retList)
}
