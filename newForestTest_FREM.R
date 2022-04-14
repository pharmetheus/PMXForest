library(dplyr)
library(ggplot2)
library(ggpubr)
library(PMXFrem)


set.seed(123)
runno   <- "22-3"
modDevDir <- "/PMX/Projects/Pharmetheus/PMX-REP-PMX-2/Analysis/Model/SimVal/"

extFile <- paste0(modDevDir,"run",runno,".ext")
covFile <- paste0(modDevDir,"run",runno,".cov")

bsFile <- paste0(modDevDir,"bs",runno,".dir/raw_results_run",runno,"bs.csv")

dfSamplesCOVfrem     <- getSamples(covFile,extFile=extFile,n=175)
dfSamplesBSfrem     <- getSamples(bsFile,extFile=extFile,n=175)

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


  return(list(CL,FREL,V))
}

functionListName <- c("CL","FREL","V")
functionList<-list(paramFunction)



#22-3
dfCovs <- createInputForestData(
  list("NCIL" = c(0,1),
       "FORM" = c(0,1),
       "FOOD" = c(0,1),
       "GENO"=list("GENO1" = c(1,0,0,0),
                   "GENO3" = c(0,0,1,0),
                   "GENO4" = c(0,0,0,1)),
       "RACEL2" =c(0,1),
       "WT"   = c(64,116),
       "BMI"  = c(23,39),
       "AGE"  = c(25,61),
       "CRCL" = c(86,150),
       "AST"  = c(15,34),
       "ALT"  = c(14,43),
       "BILI"  = c(5,15),
       "SEX"  = c(1,2)),
  iMiss=-99)


modFile      <- paste0(modDevDir,"run",runno,".mod")
noBaseThetas <- 13
noCovThetas  <-  11



covnames  <- c("NCI=0","NCI>0","Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian",
               "Other","WT 64 kg","WT 116 kg",
               "BMI 23 kg/m^2","BMI 39 kg/m^2",
               "Age 25 y","Age 61 y","CRCL 83 mL/min","CRCL 150 mL/min","AST 15 IU","AST 34 IU","ALT 14 IU","ALT 43 IU",
               "Bilirubin 5 µmol/L", "Bilirubin 15 µmol/L","Male","Female")

covariates <- c("NCI","Formulation","Food status","2D6 genotype","Race","Weight","BMI","Age","Createnine\nclearance","AST","ALT","Bilirubin","Sex")

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

forestPlot(dfresCOVfrem,plotRelative=FALSE)

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

forestPlot(dfresBSfrem,plotRelative=FALSE)
