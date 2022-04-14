library(dplyr)
library(ggplot2)
library(ggpubr)
## This file sets up the reference objects for the ubuts tests
set.seed(865765)

dfCovs <- createInputForestData(
  list("NCIL" = c(0,1),
       "FORM" = c(0,1),
       "FOOD" = c(0,1),
       "GENO"=list("GENO1" = c(1,0,0,0),
                   "GENO3" = c(0,0,1,0),
                   "GENO4" = c(0,0,0,1)),
       "RACEL"= c(1,2,3),
       "WT"   = c(65,115),
       "AGE"  = c(27,62),
       "CRCL" = c(83,150),
       "SEX"  = c(1,2)),
  iMiss=-99)

covnames  <- c("NCI=0","NCI>0","Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian",
               "African American","Asian and other","WT 65 kg","WT 115 kg",
               "Age 27 y","Age 62 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

## The print names of the covariates
covariates <- c("NCI","Formulation","Food status","2D6 genotype","Race","Weight","Age","Createnine\nclearance","Sex")

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


runno   <- 7
extFile <- system.file("extdata",paste0("newExample/run",runno,".ext"),package="PMXForest")
covFile <- system.file("extdata",paste0("newExample/run",runno,".cov"),package="PMXForest")
bsFile  <- system.file("extdata",paste0("newExample/bs",runno,".dir/raw_results_run",runno,"bs.csv"),package="PMXForest")

dfSamplesCOVscm <- getSamples(covFile,extFile=extFile,n=175)
dfSamplesBSscm  <- getSamples(bsFile,extFile=extFile,n=175)

noBaseThetas <- 14

dfresCOVscm <- getForestDFSCM(dfCovs           = dfCovs,
                              cdfCovsNames     = covnames,
                              functionList     = list(paramFunction),
                              functionListName = functionListName,
                              noBaseThetas     = noBaseThetas,
                              dfParameters     = dfSamplesCOVscm,
                              dfRefRow         = NULL
)

forestPlot(dfresCOVscm,plotRelative=FALSE,parameters = c("CL","Frel","V"))

dfresBSscm <- getForestDFSCM(dfCovs           = dfCovs,
                              cdfCovsNames     = covnames,
                              functionList     = list(paramFunction),
                              functionListName = functionListName,
                              noBaseThetas     = noBaseThetas,
                              dfParameters     = dfSamplesBSscm,
                              dfRefRow         = NULL
)

forestPlot(dfresBSscm,plotRelative=FALSE,parameters = c("CL","Frel","V"))
