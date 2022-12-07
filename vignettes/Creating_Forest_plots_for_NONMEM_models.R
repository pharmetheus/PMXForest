## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
packageName <- "PMXForest"
suppressPackageStartupMessages(library(packageName,character.only = TRUE))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

theme_set(theme_bw(base_size=18))
theme_update(plot.title = element_text(hjust = 0.5))



set.seed(865765)

## -----------------------------------------------------------------------------
dfCovs <- createInputForestData(
  list("NCIL" = c(0,1),
       "FORM" = c(0,1),
       "FOOD" = c(0,1),
       "GENO"=list("GENO1" = c(1,0,0,0),
                   "GENO3" = c(0,0,1,0),
                   "GENO4" = c(0,0,0,1)),
       "RACEL2"= c(0,1),
       "WT"   = c(65,115),
       "AGE"  = c(27,62),
       "CRCL" = c(83,150),
       "SEX"  = c(1,2)),
  iMiss=-99)

dfCovs

## -----------------------------------------------------------------------------
covnames  <- c("NCI=0","NCI>0","Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM",
               "Caucasian","Other","WT 65 kg","WT 115 kg",
               "Age 27 y","Age 62 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

## -----------------------------------------------------------------------------
covariateGroupNames <- c("NCI","Formulation","Food status","2D6 genotype","Race","Weight","Age",
                         "Createnine\nclearance","Sex")

## -----------------------------------------------------------------------------
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

  FREL <- TVFREL
  CL   <- TVCL
  V    <- TVV
  
  ## Compute AUC after 80 mg dose
  AUC <- 80/(CL/FREL)

  return(list(CL,FREL,AUC,V))
}

functionListName <- c("CL","Frel","AUC","V")

## -----------------------------------------------------------------------------
runno   <- "7"
modDir <- "SimVal"
extFile <- system.file("extdata",paste0(modDir,"/run",runno,".ext"),package=packageName)
covFile <- system.file("extdata",paste0(modDir,"/run",runno,".cov"),package=packageName)
bsFile  <- system.file("extdata",paste0(modDir,"/bs",runno,".dir/raw_results_run",runno,"bs.csv"),package=packageName)
sirFile <- system.file("extdata",paste0(modDir,"/sir",runno,".dir/raw_results_run",runno,".csv"),package=packageName)

## -----------------------------------------------------------------------------
#Covariance matrix
dfSamplesCOV <- getSamples(covFile,extFile,n=175)

#Bootstrap 
dfSamplesBS <- getSamples(bsFile,extFile)

#SIRFile
dfSamplesSIR <- getSamples(sirFile,extFile)

#Small bootstrap mvnorm sampling
dfSamplesBSsmall <- getSamples(bsFile,extFile,n=175)

## -----------------------------------------------------------------------------
dfresCOV <- getForestDFSCM(dfCovs           = dfCovs,
                           cdfCovsNames     = covnames,
                           functionList     = list(paramFunction),
                           functionListName = functionListName,
                           noBaseThetas     = 14,
                           dfParameters     = dfSamplesCOV
)
head(dfresCOV)

## -----------------------------------------------------------------------------
dfresBS <- getForestDFSCM(dfCovs           = dfCovs,
                          cdfCovsNames     = covnames,
                          functionList     = list(paramFunction),
                          functionListName = functionListName,
                          noBaseThetas     = 14,
                          dfParameters     = dfSamplesBS
)

## ----fig.width=14,fig.height=12-----------------------------------------------
forestPlot(dfresCOV,parameters="CL",groupNameLabels = covariateGroupNames)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
dataFile <- system.file("extdata",paste0(modDir,"/DAT-1-MI-PMX-2.csv"),package=packageName)
dfData   <- read.csv(dataFile) %>% distinct(ID,.keep_all=TRUE) %>% slice(1:300) # Will only use 300 to save time

## ----cache=TRUE---------------------------------------------------------------
dfresCOVemp <- getForestDFemp( dfData             = dfData,
                               covExpressionsList = lsExpr,
                               cdfCovsNames       = covnamesEmp,
                               functionList       = list(paramFunction),
                               functionListName   = functionListName,
                               noBaseThetas       = 14,
                               dfParameters       = dfSamplesCOV,
                               ncores             = 1,
                               cstrPackages       = "dplyr"
)


## ----fig.width=14,fig.height=12-----------------------------------------------
forestPlot(dfresCOVemp,parameters="CL",groupNameLabels = covariateGroupNames,noVar = FALSE)

