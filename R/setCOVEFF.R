#' setCOVEFF
#'
#' @description Set the COEFF column in dfes to TRUE based on user input. Used to specify which covariates that are significant in an SCM.
#' @inheritParams setupForestPlotData
#'
#' @param lsTrue List of two element vectors that specify parameter and covariate combination to set to TRUE.
#'
#' @return A data frame similar to \code{dfres} except that the COVEFF column has been modified.
#'
#' @examples
#' \dontrun{
#' lsTrue <- list(
#'   c("CL","FOOD"),
#'   c("CL","GENO"),
#'   c("CL","WT"),
#'   c("Frel","FORM"),
#'   c("Frel","GENO"),
#'   c("Frel","SEX"),
#'   c("AUC","FORM"),
#'   c("AUC","FOOD"),
#'   c("AUC","GENO"),
#'   c("AUC","WT"),
#'   c("AUC","SEX"),
#'   c("V","WT")
#' )
#' newDfres <- setCOVEFF(dfresCOVscm,lsTrue)
#' }
setCOVEFF <- function(dfres,lsTrue) {

  if(!class(lsTrue)=="list") strop("lsTrue must be a list")

  ## Set all COVEFF to FALSE before specifying the true ones
  dfres$COVEFF <- FALSE

  for(myVec in lsTrue)   {
    dfres <- dfres %>% mutate(COVEFF = ifelse(PARAMETER==myVec[1] & GROUPNAME== myVec[2],TRUE,COVEFF))
  }

  return(dfres)
}


