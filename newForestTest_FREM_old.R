library(dplyr)
library(ggplot2)
library(ggpubr)
library(PMXFrem)


set.seed(123)
runno   <- 62
extFile <- system.file("extdata",paste0("run",runno,".ext"),package="PMXForest")
covFile <- system.file("extdata",paste0("run",runno,".cov"),package="PMXForest")

dfSamplesCOVfrem     <- getSamples(covFile,extFile=extFile,n=200)

paramFunction <- function(basethetas, covthetas, dfrow, ...) {

  FRELGENO <- 1
  if (dfrow$GENO!=-99 && dfrow$GENO == 2) FRELGENO <- 1
  if (dfrow$GENO!=-99 && dfrow$GENO == 3) FRELGENO <- (1 + basethetas[14])
  if (dfrow$GENO!=-99 && dfrow$GENO == 4) FRELGENO <- (1 + basethetas[15])
  if (dfrow$GENO!=-99 && dfrow$GENO == 1) FRELGENO <- (1 + basethetas[16])

  FRELFORM <- 1
  if(any(names(dfrow)=="FORM") && dfrow$FORM!=-99) {
    if (dfrow$FORM == 1) FRELFORM <- 1
    if (dfrow$FORM == 0) FRELFORM <- (1 + basethetas[13])
  }

  FRELCOVTIME <- FRELFORM
  FRELCOV     <- FRELGENO

  CLFOOD <- 1
  if (dfrow$FOOD!=-99 && dfrow$FOOD == 1) CLFOOD <- 1
  if (dfrow$FOOD!=-99 && dfrow$FOOD == 0) CLFOOD <- (1 + basethetas[12])

  CLWT <- 1
  if (any(names(dfrow)=="WT") && dfrow$WT!=-99)  CLWT <- (dfrow$WT / 75)**basethetas[2]

  CLGENO <- 1
  if (dfrow$GENO!=-99 && dfrow$GENO == 2) CLGENO <- 1
  if (dfrow$GENO!=-99 && dfrow$GENO == 1) CLGENO <- (1 + basethetas[9])
  if (dfrow$GENO!=-99 && dfrow$GENO == 3) CLGENO <- (1 + basethetas[10])
  if (dfrow$GENO!=-99 && dfrow$GENO == 4) CLGENO <- (1 + basethetas[11])

  CLCOVTIME <- CLFOOD
  CLCOV     <- CLWT * CLGENO

  VWT <- 1
  if (any(names(dfrow)=="WT") && dfrow$WT!=-99) VWT <- (dfrow$WT/75)**basethetas[3]
  VCOV <- VWT

  TVFREL <- basethetas[1] * FRELCOV
  TVCL   <- basethetas[4] * CLCOV
  TVV    <- basethetas[5] * VCOV
  TVMAT  <- basethetas[6]
  TVD1   <- basethetas[7]
  TVRUV  <- basethetas[8]

  MU_1 <- log(TVRUV)
  MU_2 <- TVD1
  MU_3 <- log(TVFREL)
  MU_4 <- log(TVCL)
  MU_5 <- log(TVV)
  MU_6 <- log(TVMAT)

  RUV   <- exp(MU_1)
  D1FR  <- MU_2
  FREL  <- FRELCOVTIME * exp(MU_3 + covthetas[1])
  CL    <- CLCOVTIME   * exp(MU_4 + covthetas[2])
  V     <- exp(MU_5               + covthetas[3])
  MAT   <- exp(MU_6               + covthetas[4])

  D1 <- MAT * (1 - D1FR)
  F1 <- FREL
  KA <- 1 / (MAT - D1)


  return(list(CL,FREL,V))
}

functionListName <- c("CL","FREL","V")
functionList<-list(paramFunction)

dfCovs <- createInputForestData(
  list("FORM" = c(0,1),
       "FOOD" = c(0,1),
       "GENO" = c(1,2,3,4),
       "RACEL2" =c(0,1),
       "WT"   = c(50,75,100),
       "BMI"  = c(20,25,30),
       "AGE"  = c(50,70,85),
       "CRCL" = c(83,150),
       "SEX"  = c(1,2)),
  iMiss=-99)

modFile      <- system.file("extdata",paste0("run",runno,".mod"),package="PMXForest")
noBaseThetas <- 16
noCovThetas  <-  5

covnames  <- c("Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian",
               "Other","WT 50 kg","WT 75 kg","WT 100 kg",
               "BMI 20 kg/m^2","BMI 25 kg/m^2","BMI 30 kg/m^2",
               "Age 50 y","Age 70 y","Age 85 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

covariates <- c("Formulation","Food status","2D6 genotype","Race","Weight","BMI","Age","Createnine\nclearance","Sex")

dfresCOVfrem <-getForestDFFREM(dfCovs = dfCovs,
                               cdfCovsNames = covnames,
                               covNames=getCovNames(modFile),
                               functionList = list(paramFunction),
                               functionListName = functionListName,
                               noBaseThetas = noBaseThetas,
                               noParCov = 4,
                               noSigmas = 2,
                               noSkipOm = 2,
                               noCovThetas = noCovThetas,
                               dfParameters = dfSamplesCOVfrem,
                               probs = c(0.05, 0.95),
                               dfRefRow = NULL,
                               quiet = TRUE,
                               ncores = 1,
                               cstrPackages = c("PMXFrem","dplyr"))

forestPlot(dfresCOVfrem,groupNameLabels = covariates,plotRelative=FALSE)


#################################

set.seed(123)
runno   <- 61
extFile <- system.file("extdata",paste0("run",runno,".ext"),package="PMXForest")
bootFile<- system.file("extdata",paste0("bs",runno,".dir/raw_results_run",runno,"bs.csv"),package="PMXForest")

dfSamplesBSfrem      <- getSamples(bootFile,extFile=extFile)

paramFunction <- function(basethetas, covthetas, dfrow, ...) {

  FRELGENO <- 1
  if (dfrow$GENO!=-99 && dfrow$GENO == 2) FRELGENO <- 1
  if (dfrow$GENO!=-99 && dfrow$GENO == 3) FRELGENO <- (1 + basethetas[14])
  if (dfrow$GENO!=-99 && dfrow$GENO == 4) FRELGENO <- (1 + basethetas[15])
  if (dfrow$GENO!=-99 && dfrow$GENO == 1) FRELGENO <- (1 + basethetas[16])

  FRELFORM <- 1
  if(any(names(dfrow)=="FORM") && dfrow$FORM!=-99) {
    if (dfrow$FORM == 1) FRELFORM <- 1
    if (dfrow$FORM == 0) FRELFORM <- (1 + basethetas[13])
  }

  FRELCOVTIME <- FRELFORM
  FRELCOV     <- FRELGENO

  CLFOOD <- 1
  if (dfrow$FOOD!=-99 && dfrow$FOOD == 1) CLFOOD <- 1
  if (dfrow$FOOD!=-99 && dfrow$FOOD == 0) CLFOOD <- (1 + basethetas[12])

  CLWT <- 1
  if (any(names(dfrow)=="WT") && dfrow$WT!=-99)  CLWT <- (dfrow$WT / 75)**basethetas[2]

  CLGENO <- 1
  if (dfrow$GENO!=-99 && dfrow$GENO == 2) CLGENO <- 1
  if (dfrow$GENO!=-99 && dfrow$GENO == 1) CLGENO <- (1 + basethetas[9])
  if (dfrow$GENO!=-99 && dfrow$GENO == 3) CLGENO <- (1 + basethetas[10])
  if (dfrow$GENO!=-99 && dfrow$GENO == 4) CLGENO <- (1 + basethetas[11])

  CLCOVTIME <- CLFOOD
  CLCOV     <- CLWT * CLGENO

  VWT <- 1
  if (any(names(dfrow)=="WT") && dfrow$WT!=-99) VWT <- (dfrow$WT/75)**basethetas[3]
  VCOV <- VWT

  TVFREL <- basethetas[1] * FRELCOV
  TVCL   <- basethetas[4] * CLCOV
  TVV    <- basethetas[5] * VCOV
  TVMAT  <- basethetas[6]
  TVD1   <- basethetas[7]
  TVRUV  <- basethetas[8]

  MU_1 <- log(TVRUV)
  MU_2 <- TVD1
  MU_3 <- log(TVFREL)
  MU_4 <- log(TVCL)
  MU_5 <- log(TVV)
  MU_6 <- log(TVMAT)

  RUV   <- exp(MU_1)
  D1FR  <- MU_2
  FREL  <- FRELCOVTIME * exp(MU_3 + covthetas[1])
  CL    <- CLCOVTIME   * exp(MU_4 + covthetas[2])
  V     <- exp(MU_5               + covthetas[3])
  MAT   <- exp(MU_6               + covthetas[4])

  D1 <- MAT * (1 - D1FR)
  F1 <- FREL
  KA <- 1 / (MAT - D1)


  return(CL)
}

functionListName <- c("CL")
functionList<-list(paramFunction)

dfCovs <- createInputForestData(
  list("FORM" = c(0,1),
       "FOOD" = c(0,1),
       "GENO" = c(1,2,3,4),
       "RACE"=list("RACEL_3"=c(1,0,0),
                   "RACEL_2"=c(0,1,0)),
       "WT"   = c(50,75,100),
       "BMI"  = c(20,25,30),
       "AGE"  = c(50,70,85),
       "CRCL" = c(83,150),
       "SEX"  = c(1,2)),
  iMiss=-99)

modFile      <- system.file("extdata",paste0("run",runno,".mod"),package="PMXForest")
noBaseThetas <- 16
noCovThetas  <-  6


covnames  <- c("Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian",
               "African American","Asian and other","WT 50 kg","WT 75 kg","WT 100 kg",
               "BMI 20 kg/m^2","BMI 25 kg/m^2","BMI 30 kg/m^2",
               "Age 50 y","Age 70 y","Age 85 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

## The print names of the covariates
covariates <- c("Formulation","Food status","2D6 genotype","Race","Weight","BMI","Age","Createnine\nclearance","Sex")

dfresCOVfremBS <-getForestDFFREM(dfCovs = dfCovs,
                                            cdfCovsNames = covnames,
                                            covNames=getCovNames(modFile),
                                            functionList = list(paramFunction),
                                            functionListName = functionListName,
                                            noBaseThetas = noBaseThetas,
                                            noParCov = 4,
                                            noSigmas = 2,
                                            noSkipOm = 2,
                                            noCovThetas = noCovThetas,
                                            dfParameters = dfSamplesBSfrem ,
                                            probs = c(0.05, 0.95),
                                            dfRefRow = NULL,
                                            quiet = TRUE,
                                            ncores=6,
                                            cstrPackages = c("PMXFrem","dplyr"))


forestPlot(dfresCOVfremBS,groupNameLabels = covariates )


