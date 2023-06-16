#' createInputForestData
#'
#' @param listCovs a list with covariate names and values of the form; list("AGE"=c(1,2)) etc.
#' For multiple covariates (i.e. not univariate) a list of lists should be used, e.g. list("AGEBW"=list("AGE"=c(1,2),"BW"=c(3,4)))
#' @param iMiss The value that should be filled for missing, default = -99
#'
#' @return a data frame with a covariate value (univariate) per row or several covariate values per row (multivariate). The data frame also contains a COVARIATEGROUPS column based on the names in the input list.
#' @export
#' @examples
#' \dontrun{
#' ## One covariate effect per 'group' on the Forest plot
#' dfCovs <- createInputForestData(
#' list("FORM" = c(0,1),
#'      "FOOD" = c(0,1),
#'      "GENO" = c(1,2,3,4),
#'      "RACEL"= c(1,2,3),
#'      "WT"   = c(65,115),
#'      "AGE"  = c(27,62),
#'      "CRCL" = c(83,150),
#'      "SEX"  = c(1,2)
#'      ),
#'      Miss=-99)
#' }
#'
#' ## In this case the GENO covariate levels (2,3 and 4) will be grouped in the Forest plot.
#' dfCovs <- createInputForestData(
#' list("FORM" = c(0,1),
#'      "FOOD" = c(0,1),
#'      "GENO"=list("GENO_2" = c(0,1,0,0),
#'                  "GENO_3" = c(0,0,1,0),
#'                  "GENO_4" = c(0,0,0,1))
#'           ),
#'           iMiss=-99
#' )
#'
createInputForestData <- function(listCovs, iMiss = -99) {
  fixedname<-"COVARIATEGROUPS"
  cstrUniqueCols <- c()

    for (i in 1:length(listCovs)) {
      if (is.list(listCovs[[i]])) {
        for (j in 1:length(listCovs[[i]])) {
          cstrUniqueCols <- c(cstrUniqueCols, names(listCovs[[i]][j]))
        }
      } else {
        cstrUniqueCols <- c(cstrUniqueCols, names(listCovs[i]))
      }
    }


  cstrUniqueCols <- unique(cstrUniqueCols)
  ncols <- length(cstrUniqueCols)
  df <- data.frame(rbind(1:(ncols+1)))
  names(df) <- c(cstrUniqueCols,fixedname)
  r <- 1
  for (i in 1:length(listCovs)) {
    if (is.list(listCovs[[i]])) {
      strNameBig<-names(listCovs[i])
      strVal1 <- listCovs[[i]][[1]]
      for (k in 1:length(strVal1)) {
        dfr <- df[1, ]
        dfr[, ] <- iMiss
        for (j in 1:length(listCovs[[i]])) {
          strName <- names(listCovs[[i]][j])
          if (!is.na(listCovs[[i]][[j]][k])) dfr[[strName]] <- listCovs[[i]][[j]][k]
        }
        dfr[[fixedname]]<-strNameBig
        df <- rbind(df, dfr)
      }
    } else {
      vecr <- listCovs[[i]]
      strName <- names(listCovs[i])
      for (k in 1:length(vecr)) {
        dfr <- df[1, ]
        dfr[, ] <- iMiss
        dfr[[strName]] <- vecr[k]
        dfr[[fixedname]]<-strName
        df <- rbind(df, dfr)
      }
    }
  }
  df<-df[-1,]
  return(df)
}
