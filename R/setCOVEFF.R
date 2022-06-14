setCOVEFF <- function(dfres,lsTrue) {

  if(!class(lsTrue)=="list") strop("lsTrue must be a list")

  ## Set all COVEFF to FALSE before specifying the true ones
  dfres$COVEFF <- FALSE

  for(myVec in lsTrue)   {
    dfres <- dfres %>% mutate(COVEFF = ifelse(PARAMETER==myVec[1] & GROUPNAME== myVec[2],TRUE,COVEFF))
  }

  return(dfres)
}


# lsTrue <- list(
#   c("CL","FOOD"),
#   c("CL","GENO"),
#   c("CL","WT"),
#   c("Frel","FORM"),
#   c("Frel","GENO"),
#   c("Frel","SEX"),
#   c("AUC","FORM"),
#   c("AUC","FOOD"),
#   c("AUC","GENO"),
#   c("AUC","WT"),
#   c("AUC","SEX"),
#   c("V","WT")
#
# )
#
# nrow(setCOVEFF(dfresEMP,list()) %>% filter(COVEFF) %>% distinct(PARAMETER,GROUPNAME) %>% arrange(PARAMETER))
# setCOVEFF(dfresCOVscm,lsTrue) %>% filter(COVEFF) %>% distinct(PARAMETER,GROUPNAME) %>% arrange(PARAMETER)
