#' getCovStats
#'
#' @description To quickly get a dfCovs data from for Forest plots by just supplying a data file and a vector of covariate names.
#' @param data The data file that includes covariates.
#' @param covariates A character vector of covariate names as they appear in the data file.
#' @param minLevels The maximum number of unique values for a covariate to be regarded as continus. Devfault is 10.
#' @param probs A vector of two number - the lower and upper percentiles to use for summarising continuous covariates.
#' @param idVar The name of the ID column to use for identifying one record per subject.
#' @param missVal The missing value indicator. Default is -99.
#'
#' @return A data frame with covariate statistics to be used as an argument for Forest plots.
#' @export
#'
#' @examples
#' \dontrun{
#' modDevDir <- "./"
#' fremRun   <- 31
#' fremFiles <- PMXFrem::getFileNames(31,modDevDir = modDevDir)
#' covNames  <- PMXFrem::getCovNames(fremFiles$mod)
#'
#' data <- fread(file.path(modDevDir,"DAT-2-MI-PMX-2-onlyTYPE2-new.csv"))
#'
#' dfCovs <- createInputForestData(
#'   getCovStats(data,covNames$orgCovNames,missVal=-99),
#'   iMiss    =-99
#' )}
getCovStats <- function(data,covariates,minLevels = 10,probs=c(0.05,0.95),idVar = "ID",missVal=-99,nsig=3) {

  data <- data %>% distinct(!!ensym(idVar),.keep_all = TRUE)

  ## Remove duplicated idVar rows
  data <- subset(data,!duplicated(idVar))

  retList <- list()
  for(myCov in covariates) {
    data <- data[data[[myCov]] != missVal,]

    if(length(unique(data[[myCov]])) <= minLevels) {

      numLevs <- length(unique(data[[myCov]]))
      ## Binary covariate
      if(numLevs == 2) {
        retList[[myCov]] <- unique(data[[myCov]])
      } else {

        levs <- sort(unique(data[[myCov]]))
        covList <- list()
        for (i in 2:numLevs) {

          vec <- rep(0,numLevs)
          vec[i] <- 1

          covList[[paste0(myCov,"_",levs[i])]] <- vec
        }

        retList[[myCov]] <- covList

      }

    } else {
      retList[[myCov]] <- round(quantile(data[[myCov]],p=probs),digits = nsig)
    }
  }
  return(retList)
}
