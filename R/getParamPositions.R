#' getParamPositions
#'
#' @description Figure out the parameter positions in a PsN raw_results file.
#'
#' @importFrom stringr str_match
#' @param fileName The name of a PsN raw_results_structure file.
#'
#' @return A list of parameter positions THETA for thetas, SIGMA for sigmas and OMEGA for omegas.
#' @export
#'
getParamPositions <- function(fileName) {

  con <- file(fileName, "r")
  thetaid <- NULL
  omegaid <- NULL
  sigmaid <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    omstr <- str_match(line, "^omega=([0-9]+),([0-9]+)")
    if (!is.na(omstr[, 1])) {
      if (!is.null(omegaid)) break
      omegaid <- as.numeric(omstr[, 2])
      omeganum <- as.numeric(omstr[, 3])
    }
    thstr <- str_match(line, "^theta=([0-9]+),([0-9]+)")
    if (!is.na(thstr[, 1])) {
      if (!is.null(thetaid)) break
      thetaid <- as.numeric(thstr[, 2])
      thetanum <- as.numeric(thstr[, 3])
    }
    sigstr <- str_match(line, "^sigma=([0-9]+),([0-9]+)")
    if (!is.na(sigstr[, 1])) {
      if (!is.null(sigmaid)) break
      sigmaid <- as.numeric(sigstr[, 2])
      sigmanum <- as.numeric(sigstr[, 3])
    }

    if (!is.null(omegaid) & !is.null(thetaid) & !is.null(sigmaid)) {
      break
    }
  }

  close(con)

  tmp <- NULL
  if (!is.null(thetaid)) tmp$THETA <- c(thetaid, thetanum)
  if (!is.null(sigmaid)) tmp$SIGMA <- c(sigmaid, sigmanum)
  if (!is.null(omegaid)) tmp$OMEGA <- c(omegaid, omeganum)
  return(tmp)
}
