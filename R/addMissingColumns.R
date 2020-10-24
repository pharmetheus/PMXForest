#' addMissingColumns
#'
#' @description Transforms a PsN raw reasults file to a format that is the same as a NONMEM ext file. This involves adding
#' columns for fixed parameters, rearranging and renaming columns.
#'
#' @param dfParams The file with a PsN raw results structure.
#' @param dfExt The NONMEM ext file that is used as a temlate for the rearranging
#' @param zerosindex A vector with indicies that indicates which parameters that are fixed to zero.
#'
#' @return A data frame with a structure similar to a NONMEM ext file.
#' @export
#'
#' @examples
#' \dontrun{
#' dfParameters <- read.csv(file = input, header = TRUE)
#' dfExt <- subset(getExt(extFile = extFile), ITERATION == "-1000000000")
#' dfParameters <- addMissingColumns(dfParameters, dfExt)
#' }
#'
addMissingColumns <- function(dfParams, dfExt, zerosindex) {
  if (is.null(zerosindex)) {
    # Get the OMGEA/SIGMA off-diag parameters which are zero and hence assume they are not included in dfParams
    # zerosindex<-which(grepl("^(OMEGA|SIGMA)\\.{1}(?!.*([0-9]+)\\.{1}\\2).*\\.{1}",names(dfExt),perl=T) & dfExt[1,]==0)
    zerosindex <- NULL
    tmpz <- (dfExt[1, ] == 0)
    namesdf <- names(dfExt)
    for (i in 1:length(dfExt)) {
      for (j in 1:i) {
        if (i != j) zerosindex <- c(zerosindex, which(tmpz & grepl(paste0("^(OMEGA|SIGMA)\\.", i, "\\.", j, "\\."), namesdf)))
      }
    }
  }

  if (!is.null(zerosindex) && length(zerosindex)>0) {
    # Get the parameters which are zero
    dft <- dfExt[rep(1, nrow(dfParams)), zerosindex] # Repeat to the same size as dfParams
    names(dfParams) <- names(dfExt)[-c(1, zerosindex, ncol(dfExt))]

    # Add the cols to dfParams
    dfParams <- cbind(dfParams, dft)

    # Sort dfParams correctly
    dfParams <- dfParams[, names(dfExt)[-c(1, ncol(dfExt))]]
  } else {
    names(dfParams)<-names(dfExt)[-c(1,ncol(dfExt))]
  }
  return(dfParams)
}
