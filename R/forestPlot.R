#' Get the correct dfres variables to plot based on user specification
#'
#' @description Determine which variables from dfres to base the forest plot on
#' @keywords internal
#' @param plotRelative Should the plot be made on the relative scale (TRUE or FALSE)
#' @param noVar Should the uncertainty in the reference be included (FALSE) or not (TRUE)
#' @param reference Which reference should be used, the final parameter estimates (final) or the parameter function (func)
#'
#' @return A vector with the names from dfres to use for POINT, q1, q2 and ref
#' @export
#'
#' @examples
#' getPlotVars()
#'
getPlotVars <- function(plotRelative=TRUE,noVar=FALSE,reference="func") {

  if(!(reference %in% c("func","final"))) {stop("reference needs to be either func or final.")}

  if(!plotRelative & noVar & reference == "func") {
    vars <- c(point="POINT",q1="Q1",q2="Q2",ref="REFFUNC")

  } else if(!plotRelative & noVar & reference == "final") {
    vars <- c(point="POINT",q1="Q1",q2="Q2",ref="REFTRUE")

  } else if(!plotRelative & !noVar & reference == "func") {
    stop("The combination of plotRelative=FALSE and noVar=TRUE is not possible\n.")

  } else if(!plotRelative & !noVar & reference == "final") {
    stop("The combination of plotRelative=FALSE and noVar=TRUE is not possible\n.")

  } else if(plotRelative & noVar & reference == "func") {
    vars <- c(point="POINT_NOVAR_REL_REFFUNC",q1="Q1_NOVAR_REL_REFFUNC",q2="Q2_NOVAR_REL_REFFUNC",ref="REFFUNC")

  } else if(plotRelative & noVar & reference == "final") {
    vars <- c(point="POINT_NOVAR_REL_REFTRUE",q1="Q1_NOVAR_REL_REFTRUE",q2="Q2_NOVAR_REL_REFTRUE",ref="REFTRUE")

  } else if(plotRelative & !noVar & reference == "func") {
    vars <- c(point="POINT_REL_REFFUNC",q1="Q1_REL_REFFUNC",q2="Q2_REL_REFFUNC",ref="REFFUNC")

  } else if(plotRelative & !noVar & reference == "final") {
    vars <- c(point="POINT_REL_REFTRUE",q1="Q1_REL_REFTRUE",q2="Q2_REL_REFTRUE",ref="REFTRUE")

  }
  return(vars)
}
