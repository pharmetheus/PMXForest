library(dplyr)
library(ggplot2)
library(ggpubr)


set.seed(865765)

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
modDir <- "SimVal"
extFile <- system.file("extdata",paste0(modDir,"/run",runno,".ext"),package="PMXForest")
covFile <- system.file("extdata",paste0(modDir,"/run",runno,".cov"),package="PMXForest")

dfSamplesCOVscm <- getSamples(covFile,extFile=extFile,n=175)

noBaseThetas <- 14

######################################
## Calculate Forest plot statistics ##
######################################
dfresCOVscm <- getForestDFSCM(dfCovs           = dfCovs,
                              cdfCovsNames     = covnames,
                              functionList     = list(paramFunction),
                              functionListName = functionListName,
                              noBaseThetas     = noBaseThetas,
                              dfParameters     = dfSamplesCOVscm,
                              dfRefRow         = NULL
)

#####################
## Create the plot ##
#####################

forestPlot(dfresCOVscm,plotRelative=TRUE,noVar=FALSE,parameters = c("CL"),sigdigits = 3)

forestPlot(dfresCOVscm,plotRelative=TRUE,noVar=FALSE,parameters = c("CL"),sigdigits = 3,onlySignificant = TRUE)

forestPlot(dfresCOVscm,plotRelative=TRUE,noVar=FALSE,parameters = c("CL","Frel"),sigdigits = 3,stackedPlots = TRUE,keepYlabs = TRUE,rightStrip = FALSE,keepRightStrip = TRUE,table = FALSE)

#############################################################
## Create the same plot based on the bootstrap information ##
#############################################################

bsFile         <- system.file("extdata",paste0(modDir,"/bs",runno,".dir/raw_results_run",runno,"bs.csv"),package="PMXForest")
dfSamplesBSscm <- getSamples(bsFile,extFile=extFile,n=175)

dfresBSscm <- getForestDFSCM(dfCovs           = dfCovs,
                              cdfCovsNames     = covnames,
                              functionList     = list(paramFunction),
                              functionListName = functionListName,
                              noBaseThetas     = noBaseThetas,
                              dfParameters     = dfSamplesBSscm,
                              dfRefRow         = NULL
)

forestPlot(dfresBSscm,plotRelative=TRUE,noVar=FALSE,parameters = c("CL"),sigdigits = 3)
