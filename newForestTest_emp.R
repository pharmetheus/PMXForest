library(dplyr)
library(ggplot2)
library(ggpubr)

packageName <- "PMXForest"
set.seed(865765)

################################
## Define the covariate input ##
################################

lsExpr<-list("NCI" = expression(NCIL==0),
             "NCI" = expression(NCIL==1),
             "FORM"= expression(FORM==0),
             "FORM"= expression(FORM==1),
             "FOOD"= expression(FOOD==0),
             "FOOD"= expression(FOOD==1),
             "GENO"= expression(GENO==1),
             "GENO"= expression(GENO==2),
             "GENO"= expression(GENO==3),
             "GENO"= expression(GENO==4),
             "RACE"= expression(RACEL2==0),
             "RACE"= expression(RACEL2==1),
             "WT"  = expression(WT<70),
             "WT"  = expression(WT>104),
             "AGE" = expression(AGE<35),
             "AGE" = expression(AGE>=57),
             "CRCL"= expression(CRCL<94),
             "CRCL"= expression(CRCL>=146),
             "SEX" = expression(SEX==1),
             "SEX" = expression(SEX==2)
)

## The names associated with the entries in lsExpr
covnamesEmp <- c("NCI=0","NCI>0","Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM",
                 "Caucasian","Other","WT<70 kg","WT>104 kg",
                 "Age<35 y","Age>57 y","CRCL<94 mL/min","CRCL>146 mL/min","Male","Female")

## The print names of the covariates (covariate groups)
covariates <- c("NCI","Formulation","Food status","2D6 genotype","Race","Weight","Age",
                "Createnine\nclearance","Sex")

###################################
## Define the parameter function ##
###################################

paramFunction <- function(thetas, df, ...) {

  CLFOOD <- 1
  if (any(names(df)=="FOOD") && df$FOOD !=-99 && df$FOOD == 0) CLFOOD <- (1 + thetas[11])

  CLCOV <- CLFOOD

  FRELGENO4 <- 1
  if (any(names(df)=="GENO4") && df$GENO4 !=-99 && df$GENO4 == 1) FRELGENO4 <- (1 + thetas[13])

  FRELFORM <- 1
  if (any(names(df)=="FORM") && df$FORM !=-99 && df$FORM == 0) FRELFORM <- (1 + thetas[12])

  FRELSEX <- 1
  if (any(names(df)=="SEX") && df$SEX !=-99 && df$SEX == 2) FRELSEX <- (1 + thetas[14])

  FRELCOV <- FRELSEX * FRELFORM * FRELGENO4

  TVFREL <- thetas[1]

  TVFREL <- FRELCOV * TVFREL

  if(any(names(df)=="WT") && df$WT !=-99) {
    TVCL <- thetas[4] * (df$WT / 75)**thetas[2]
  } else {
    TVCL <- thetas[4]
  }

  if (any(names(df)=="GENO1") && df$GENO1 !=-99 && df$GENO1 == 1) TVCL <- TVCL * (1 + thetas[8])
  if (any(names(df)=="GENO3") && df$GENO3 !=-99 && df$GENO3 == 1) TVCL <- TVCL * (1 + thetas[9])
  if (any(names(df)=="GENO4") && df$GENO4 !=-99 && df$GENO4 == 1) TVCL <- TVCL * (1 + thetas[10])

  TVCL <- CLCOV * TVCL

  if(any(names(df)=="WT") && df$WT !=-99) {
    TVV <- thetas[5] * (df$WT / 75)**thetas[3]
  } else {
    TVV <- thetas[5]
  }

  TVMAT <- thetas[6]
  TVMAT <- TVMAT

  TVD1 <- thetas[7]

  FREL <- TVFREL
  CL   <- TVCL
  V    <- TVV
  MAT  <- TVMAT
  D1   <- MAT * (1 - TVD1)
  KA   <- 1 / (MAT - D1)

  ## Compute AUC after 80 mg dose

  AUC <- 80/(CL/FREL)

  return(list(CL,FREL,AUC,V,MAT))
}

functionListName <- c("CL","Frel","AUC","V","MAT")

####################################################
## Define the NONMEM run information and sampling ##
####################################################
runno   <- "7"
modDir  <- "SimVal"

extFile <- system.file("extdata",paste0(modDir,"/run",runno,".ext"),package=packageName)
covFile <- system.file("extdata",paste0(modDir,"/run",runno,".cov"),package=packageName)

dfSamplesCOVscm <- getSamples(covFile,extFile=extFile,n=175)

dataFile <- system.file("extdata",paste0(modDir,"/DAT-1-MI-PMX-2.csv"),package=packageName)
dfData   <- read.csv(dataFile) %>% distinct(ID,.keep_all=TRUE)

noBaseThetas <- 14

######################################
## Calculate Forest plot statistics ##
######################################
dfresEmp <- getForestDFemp(
  dfData             = dfData,
  covExpressionsList = lsExpr,
  cdfCovsNames       = covnamesEmp,
  functionList       = list(paramFunction),
  functionListName   = functionListName,
  metricFunction     = median,
  noBaseThetas       = noBaseThetas,
  dfParameters       = dfSamplesCOVscm,
  dfRefRow           = NULL,
  ncores             = 6,
  cstrPackages       = "dplyr"
)

#####################
## Create the plot ##
#####################

lsTrue <- list(
  c("CL","FOOD"),
  c("CL","GENO"),
  c("CL","WT"),
  c("Frel","FORM"),
  c("Frel","GENO"),
  c("Frel","SEX"),
  c("AUC","FORM"),
  c("AUC","FOOD"),
  c("AUC","GENO"),
  c("AUC","WT"),
  c("AUC","SEX"),
  c("V","WT")
)
forestPlot(dfresEmp,plotRelative=TRUE,noVar=FALSE,parameters = c("CL"),sigdigits = 3,onlySignificantErrorBars = TRUE,setSignEff = lsTrue)

