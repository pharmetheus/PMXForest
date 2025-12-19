#' Get the correct dfres variables to plot based on user specification
#'
#' @description Determine which variables from dfres to base the forest plot on
#'
#' @importFrom stats cov median quantile
#' @importFrom utils read.csv read.table
#' @keywords internal
#' @param plotRelative Should the plot be made on the relative scale (TRUE or FALSE)
#' @param noVar Should the uncertainty in the reference be included (FALSE) or not (TRUE)
#' @param reference Which reference should be used, the final parameter estimates (final) or the parameter function (func)
#'
#' @return A vector with the names from dfres to use for POINT, q1, q2 and ref
#'
#'
#' @examples
#' # Get the relative variables with uncertainty in the reference:
#' # "POINT_REL_REFFUNC" "Q1_REL_REFFUNC" "Q2_REL_REFFUNC" "REFFUNC"
#'
#' getPlotVars()
#'
#' # Get the relative variables without uncertainty in the reference:
#' # "POINT_NOVAR_REL_REFFUNC" "Q1_NOVAR_REL_REFFUNC" "Q2_NOVAR_REL_REFFUNC" "REFFUNC"
#'
#' getPlotVars(noVar=TRUE)
#'
getPlotVars <- function(plotRelative=TRUE,noVar=FALSE,reference="func") {

  if(!(reference %in% c("func","final"))) {stop("reference needs to be either func or final.")}

  if(!plotRelative & noVar & reference == "func") {
    vars <- c(REF="REFFUNC",point="POINT",q1="Q1",q2="Q2")

  } else if(!plotRelative & noVar & reference == "final") {
    vars <- c(REF="REFFINAL",point="POINT",q1="Q1",q2="Q2")

  } else if(!plotRelative & !noVar & reference == "func") {
    stop("The combination of plotRelative=FALSE and noVar=FALSE is not possible\n.")

  } else if(!plotRelative & !noVar & reference == "final") {
    stop("The combination of plotRelative=FALSE and noVar=FALSE is not possible\n.")

  } else if(plotRelative & noVar & reference == "func") {
    vars <- c(REF="REFFUNC",point="POINT_NOVAR_REL_REFFUNC",q1="Q1_NOVAR_REL_REFFUNC",q2="Q2_NOVAR_REL_REFFUNC")

  } else if(plotRelative & noVar & reference == "final") {
    stop("The combination of plotRelative=TRUE, noVar=TRUE and reference=final is not possible\n.")

  } else if(plotRelative & !noVar & reference == "func") {
    vars <- c(REF="REFFUNC",point="POINT_REL_REFFUNC",q1="Q1_REL_REFFUNC",q2="Q2_REL_REFFUNC")

  } else if(plotRelative & !noVar & reference == "final") {
    vars <- c(REF="REFFINAL",point="POINT_REL_REFFINAL",q1="Q1_REL_REFFINAL",q2="Q2_REL_REFFINAL")

  }
  return(vars)
}

