library(dplyr)
library(ggplot2)
library(ggpubr)
library(PMXFrem)

set.seed(865765)

modDevDir  <- "/PMX/Projects/Pharmetheus/PMX-REP-PMX-2/Analysis/Model/SimNeb-for-example"
finalRun   <- 10
finalBSRun <- 11

extFile  <- file.path(modDevDir ,paste0("run",finalRun,".ext"))
dfExt    <- getExt(extFile) %>% filter(ITERATION==-1000000000)

modFile  <- file.path(modDevDir,paste0("run",finalRun,".mod"))
covNames <- getCovNames(modFile = modFile)

covFile  <- file.path(modDevDir ,paste0("run",finalRun,".cov"))
bsFile   <- file.path(modDevDir,paste0("bs",finalBSRun,".dir/raw_results_run",finalBSRun,".csv"))

dfCovs <- createInputForestData(
  list(NCI=list("NCI_1"  = c(0,1,0),
                "NCI_2"  = c(0,0,1)),
       "FORM" = c(0,1),
       "FOOD" = c(0,1),
       GENO   = list("GENO_1" = c(1,0,0,0),
                     "GENO_3" = c(0,0,1,0),
                     "GENO_4" = c(0,0,0,1)
       ),
       RACE=list("RACEL_2" = c(0,1,0),
                 "RACEL_3" = c(0,0,1)),
       "WT"   = c(64,116),
       "BMI"  = c(23,39),
       "AGE"  = c(25,61),
       "CRCL" = c(86,150),
       "AST"  = c(15,34),
       "ALT"  = c(14,43),
       "BILI"  = c(5,15),
       "SEX"  = c(1,2)),
  iMiss=-99)

covnames <- c("NCI=0","NCI=1","NCI=2",
              "Oral tablets","FDC",
              "Fasted","Fed",
              "CYP2D6 PM","CYP2D6 IM","CYP2D6 EM","CYP2D6 UM",
              "Caucasian","African American","Asian and other",
              "WT 64 kg","WT 116 kg",
              "BMI 23 kg/m2","BMI 39 kg/m2",
              "Age 25 y","Age 61 y",
              "CRCL 86 mL/min","CRCL 150 mL/min",
              "AST 15 U/L","AST 34 U/L",
              "ALT 14 U/L","ALT 43 U/L",
              "BILI 5 umol/L","BILI 15 umol/L",
              "Male","Female")




paramFunction <- function(basethetas,covthetas, dfrow, ...) {

  #if(covthetas[1] !=0) {
  # if(dfrow$WT!=-99)
  #cat(dfrow,"\n")
  # browser()
  # }

  if(any(names(dfrow)=="WT") && dfrow$WT !=-99) {
    CL <- basethetas[4] * (dfrow$WT /75)**basethetas[1]*exp(covthetas[2])
  } else {
    CL <- basethetas[4]*exp(covthetas[2])
  }

  FRELFOOD <- 1
  if (any(names(dfrow)=="FOOD") && dfrow$FOOD !=-99 && dfrow$FOOD == 0) FRELFOOD <- (1 + basethetas[9])
  FREL <- basethetas[3]*FRELFOOD*exp(covthetas[1])

  MATFOOD <- 1
  if (any(names(dfrow)=="FOOD") && dfrow$FOOD !=-99 && dfrow$FOOD == 0) MATFOOD <- (1 + basethetas[9])
  MAT <- basethetas[6]*MATFOOD*exp(covthetas[4])

  if(any(names(dfrow)=="WT") && dfrow$WT !=-99) {
    V <- basethetas[5] * (dfrow$WT / 75)**basethetas[2]*exp(covthetas[3])
  } else {
    V <- basethetas[5]*exp(covthetas[3])
  }

  # ## Compute AUC after 80 mg dose
  AUC <- 5/(CL/FREL)

  return(list(CL,FREL,AUC,V,MAT))
}

functionListName <- c("CL","Frel","AUC","V","MAT")

dfSamplesBS <- getSamples(bsFile,extFile=extFile,n=200)

dfRefRow <- data.frame(WT=85,SEX=1,FORM=1,FOOD=1)

## Set up dfCovs to reflect the reference subject
dfCovsRef <- dfCovs %>%
  mutate(WT   = ifelse(WT==-99,85,WT)) %>%
  mutate(SEX  = ifelse(SEX==-99,1,SEX)) %>%
  mutate(FORM = ifelse(FORM==-99,1,FORM)) %>%
  mutate(FOOD = ifelse(FOOD==-99,1,FOOD))

dfres <- getForestDFFREM(
  dfCovs           = dfCovsRef,
  cGrouping        = c(1,1,1, 2,2, 3,3, 4,4,4,4, 5,5,5, 6,6, 7,7, 8,8, 9,9, 10,10, 11,11, 12,12, 13,13),
  covNames         = covNames,
  cdfCovsNames     = covnames,
  functionList     = list(paramFunction),
  functionListName = functionListName,
  noBaseThetas     = 9,
  noSkipOm         = 2,
  noParCov         = 4,
  noCovThetas      = length(covNames$covNames),
  noSigmas         = 2,
  dfParameters     = dfSamplesBS,
  cores            = 6,
  dfRefRow         = dfRefRow,
  cstrPackages     = c("dplyr", "stringr"),
  quiet            = TRUE)

forestPlot(dfres,parameters = "Frel")


#################################
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
  CLCOV     <- CLWT #* CLGENO

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


runno   <- 61
extFile <- system.file("extdata",paste0("run",runno,".ext"),package="PMXForest")
covFile <- system.file("extdata",paste0("run",runno,".cov"),package="PMXForest")
bootFile<- system.file("extdata",paste0("bs",61,".dir/raw_results_run",61,"bs.csv"),package="PMXForest")

dfSamplesBSfrem <- getSamples(bootFile,extFile=extFile,n=100)


modFile      <- system.file("extdata",paste0("run",runno,".mod"),package="PMXForest")
noBaseThetas <- 16
noCovThetas  <-  5

covnames  <- c("Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian",
               "African American","Asian and other","WT 65 kg","WT 115 kg",
               "Age 27 y","Age 62 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

## The print names of the covariates
covariates <- c("Formulation","Food status","2D6 genotype","Race","Weight","Age","Createnine\nclearance","Sex")

dfresCOVfrem <-getForestDFFREM(dfCovs = dfCovs,
                               #cdfCovsNames = covnames,
                               covNames= getCovNames(modFile),
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
                               ncores = 2,
                               cstrPackages = c("PMXFrem","dplyr"))



refRow <- dfCovs %>% mutate(
  FORM = -99,
  FOOD = -99,
  GENO = -99,
  RACEL= -99,
  WT   =- 99,
  AGE  = -99,
  CRCL = -99,
  SEX = -99) %>%
  group_by(COVARIATEGROUPS) %>%
  mutate(WT=sample(c(50,75,100),1))

dfresRefRow <- getForestDFSCM(dfCovs           = dfCovs,
                              cdfCovsNames     = covnames,
                              functionList     = list(paramFunction),
                              functionListName = functionListName,
                              noBaseThetas     = noBaseThetas,
                              dfParameters     = dfSamplesCOV,
                              dfRefRow         = refRow
)


plotData       <- setupForestPlotData(dfres,plotRelative=FALSE,noVar=TRUE)
plotDataRefRow <- setupForestPlotData(dfresRefRow,plotRelative=FALSE,noVar=TRUE,sigdig=3)


old <- theme_set(theme_bw())

forestPlot(dfres)
forestPlot(dfres,tabTextSize=5)
forestPlot(dfresRefRow,plotRelative = FALSE)
forestPlot(dfres,plotData=plotData)
forestPlot(dfres,plotRelative = FALSE)
forestPlot(dfres,parameters = c("CL","FREL"),parameterLabels=c("CL(L/h)","F[rel]"),groupNameLabels=covariates)
forestPlot(dfres,plotData=myPlotData,parameters = c("CL","FREL"),parameterLabels=c("CL(L/h)","F[rel]"),groupNameLabels=covariates)
