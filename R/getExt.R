#' getExt
#'
#' @description Extracts the NONMEM iteration information from a NONMEM .ext file.
#' @param extFile The name of the .ext file.
#' @param set The $ESTIMATION the iteration information should be extracted from. Will use the last $ESTIMATION if set to \code{NULL}.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' \dontrun{
#' ## To extract the final parameter estimates
#' dfext    <- subset(getExt(extFile = "myExtFile.ext"),ITERATION=="-1000000000")
#' }
getExt <- function(extFile,set=NULL) {

  tmp   <- scan(extFile,what="character",sep="\n",quiet=TRUE)
  tabs  <- grep("TABLE",tmp)
  if(is.null(set)) set <- length(tabs)

  if(set==1 & length(tabs)==1) { # Only one set of results
    myext <- read.table(extFile,skip=1,header=T)
  } else if(set== 1 & length(tabs)>1) {
    myext <- read.table(extFile,skip=1,nrows=tabs[2]-3,header=T)
  } else if(set==2 & length(tabs)==2) {
    myext <- read.table(extFile,skip=tabs[2],header=T)
  } else if(set==2 & length(tabs)==3) {
    myext <- read.table(extFile,skip=tabs[2],nrows=length(tmp)-tabs[2]-(length(tmp)-tabs[3])-2,header=T)
  } else if(set==3 & length(tabs)==3) {
    myext <- read.table(extFile,skip=tabs[3],header=T)
  } else if(set==4 & length(tabs)==4) {
    myext <- read.table(extFile,skip=tabs[4],header=T)
  }

  return(myext)
}
