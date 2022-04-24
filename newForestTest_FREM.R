library(dplyr)
library(ggplot2)
library(ggpubr)
library(PMXFrem)

set.seed(123)

################################
## Define the covariate input ##
################################

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

## The names associated with the entries in dfCovs
covnames  <- c("NCI=0","NCI>0","Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM",
               "Caucasian","Other","WT 65 kg","WT 115 kg",
               "Age 27 y","Age 62 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

## The print names of the covariates (covariate groups)
covariates <- c("NCI","Formulation","Food status","2D6 genotype","Race","Weight","Age",
                "Createnine\nclearance","Sex")

###################################
## Define the parameter function ##
###################################

paramFunction <- function(basethetas, covthetas, dfrow, ...) {

  FRELGENO4 <- 1
  if (dfrow$GENO4!=-99 && dfrow$GENO4 == 0) FRELGENO4 <- 1
  if (dfrow$GENO4!=-99 && dfrow$GENO4 == 1) FRELGENO4 <- (1 + basethetas[13])

  FRELFORM <- 1
  if(any(names(dfrow)=="FORM") && dfrow$FORM!=-99) {
    if (dfrow$FORM == 1) FRELFORM <- 1
    if (dfrow$FORM == 0) FRELFORM <- (1 + basethetas[12])
  }

  FRELCOVTIME <- FRELFORM
  FRELCOV     <- FRELGENO4

  CLFOOD <- 1
  if (dfrow$FOOD!=-99 && dfrow$FOOD == 1) CLFOOD <- 1
  if (dfrow$FOOD!=-99 && dfrow$FOOD == 0) CLFOOD <- (1 + basethetas[11])

  CLWT <- 1
  if (any(names(dfrow)=="WT") && dfrow$WT!=-99)  CLWT <- (dfrow$WT / 75)**basethetas[2]

  CLGENO1 <- 1
  CLGENO2 <- 1
  CLGENO3 <- 1
  CLGENO4 <- 1
  if (dfrow$GENO1!=-99 && dfrow$GENO1 == 1) CLGENO1 <- (1 + basethetas[8])
  if (dfrow$GENO3!=-99 && dfrow$GENO3 == 1) CLGENO3 <- (1 + basethetas[9])
  if (dfrow$GENO4!=-99 && dfrow$GENO4 == 1) CLGENO4 <- (1 + basethetas[10])

  CLCOVTIME <- CLFOOD
  CLCOV     <- CLWT*CLGENO1*CLGENO2*CLGENO3*CLGENO4

  VWT <- 1
  if (any(names(dfrow)=="WT") && dfrow$WT!=-99) VWT <- (dfrow$WT/75)**basethetas[3]
  VCOV <- VWT

  TVFREL <- basethetas[1] * FRELCOV
  TVCL   <- basethetas[4] * CLCOV
  TVV    <- basethetas[5] * VCOV
  TVMAT  <- basethetas[6]
  TVD1   <- basethetas[7]

  MU_2 <- TVD1
  MU_3 <- log(TVFREL)
  MU_4 <- log(TVCL)
  MU_5 <- log(TVV)
  MU_6 <- log(TVMAT)

  D1FR  <- MU_2
  FREL  <- FRELCOVTIME * exp(MU_3 + covthetas[1])
  CL    <- CLCOVTIME   * exp(MU_4 + covthetas[2])
  V     <- exp(MU_5               + covthetas[3])
  MAT   <- exp(MU_6               + covthetas[4])

  D1 <- MAT * (1 - D1FR)
  F1 <- FREL
  KA <- 1 / (MAT - D1)

  AUC <- 80/(CL/FREL)

  return(list(CL,FREL,AUC,V,MAT))
}

functionListName <- c("CL","Frel","AUC","V","MAT")


####################################################
## Define the NONMEM run information and sampling ##
####################################################
runno  <- "22-3"
modDir <- "SimVal"

extFile <- system.file("extdata",paste0(modDir,"/run",runno,".ext"),package="PMXForest")
covFile <- system.file("extdata",paste0(modDir,"/run",runno,".cov"),package="PMXForest")

dfSamplesCOVfrem  <- getSamples(covFile,extFile=extFile,n=175)

modFile      <- system.file("extdata",paste0(modDir,"/run",runno,".mod"),package="PMXForest")
noBaseThetas <- 13
noCovThetas  <- 11

######################################
## Calculate Forest plot statistics ##
######################################

dfresCOVfrem <-getForestDFFREM(dfCovs           = dfCovs,
                               cdfCovsNames     = covnames,
                               covNames         = getCovNames(modFile),
                               functionList     = list(paramFunction),
                               functionListName = functionListName,
                               noBaseThetas     = noBaseThetas,
                               noParCov         = 4,
                               noSigmas         = 2,
                               noSkipOm         = 2,
                               noCovThetas      = noCovThetas,
                               dfParameters     = dfSamplesCOVfrem,
                               probs            = c(0.05, 0.95),
                               dfRefRow         = NULL,
                               quiet            = TRUE,
                               ncores           = 1,
                               cstrPackages     = c("PMXFrem","dplyr"))

dfCovsCRCL <- createInputForestData(
  list("CRCL"   = c(30,50,80,120),
       "CRCL"  = c(84,150)),
  iMiss=-99)

dfCovsCRCL$COVARIATEGROUPS[c(5,6)] <- "CRCL2"

## The names associated with the entries in dfCovs
covnamesCRCL  <- c("CRCL 30 mL/min(*)","CRCL 50 mL/min","CRCL 80 mL/min","CRCL 120 mL/min","CRCL 84 mL/min","CRCL 150 mL/min")

## The print names of the covariates (covariate groups)
covariatesCRCL <- c("Createnine\nclearance groups","Createnine\nclearance percentiles")

dfresCOVfremCRCL <-getForestDFFREM(dfCovsCRCL,
                                   cdfCovsNames     = covnamesCRCL,
                                   covNames         = PMXFrem::getCovNames(modFile),
                                   functionList     = list(paramFunction),
                                   functionListName = functionListName,
                                   noBaseThetas     = 13,#noBaseThetas,
                                   noParCov         = 4,
                                   noSigmas         = 2,
                                   noSkipOm         = 2,
                                   noCovThetas      = 11,#noCovThetas,
                                   dfParameters     = dfSamplesCOVfrem,
                                   probs            = c(0.05, 0.95),
                                   dfRefRow         = NULL,
                                   quiet            = TRUE,
                                   ncores           = 1,
                                   cstrPackages     = c("PMXFrem","dplyr"))
#####################
## Create the plot ##
#####################

forestPlot(dfresCOVfrem,plotRelative=TRUE,noVar=FALSE,parameters = c("CL"),sigdigits = 3)


#############################################################
## Create the same plot based on the bootstrap information ##
#############################################################

bsFile          <- paste0(modDevDir,"bs",runno,".dir/raw_results_run",runno,"bs.csv")
dfSamplesBSfrem <- getSamples(bsFile,extFile=extFile,n=175)

dfresBSfrem <-getForestDFFREM(dfCovs = dfCovs,
                               cdfCovsNames = covnames,
                               covNames=getCovNames(modFile),
                               functionList = list(paramFunction),
                               functionListName = functionListName,
                               noBaseThetas = noBaseThetas,
                               noParCov = 4,
                               noSigmas = 2,
                               noSkipOm = 2,
                               noCovThetas = noCovThetas,
                               dfParameters = dfSamplesBSfrem,
                               probs = c(0.05, 0.95),
                               dfRefRow = NULL,
                               quiet = TRUE,
                               ncores = 1,
                               cstrPackages = c("PMXFrem","dplyr"))

forestPlot(dfresBSfrem,plotRelative=TRUE,noVar=FALSE,parameters = c("CL"),sigdigits = 3)
