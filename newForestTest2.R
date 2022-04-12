library(dplyr)
library(ggplot2)
library(ggpubr)
## Based on the example in the template development report

set.seed(865765)

finalRun  <- 7
modDevDir <- "/PMX/Projects/Pharmetheus/PMX-REP-PMX-2/Analysis/Model/SimVal-for-example"


dfCovs <- createInputForestData(
  list("NCIL" = c(0,1),
       "FORM" = c(0,1),
       "FOOD" = c(0,1),
       GENO=list("GENO1" = c(1,0,0,0),
                 "GENO3" = c(0,0,1,0),
                 "GENO4" = c(0,0,0,1)),
       "RACEL"= c(1,2,3),
       "WT"   = c(64,116),
       "BMI"  = c(23,39),
       "AGE"  = c(25,61),
       "CRCL" = c(86,150),
       "AST"  = c(15,34),
       "ALT"  = c(14,43),
       "BILI"  = c(5,15),
       "SEX"  = c(1,2)),
  iMiss=-99)


covnames <- c("NCI=0","NCI>0","Oral tablets","FDC","Fasted","Fed","CYP2D6 UM","CYP2D6 EM","CYP2D6 IM","CYP2D6 PM","Caucasian",
              "African American","Asian and other",
              "WT 64 kg","WT 116 kg",
              "BMI 23 kg/m2","BMI 39 kg/m2",
              "Age 25 y","Age 61 y",
              "CRCL 86 mL/min","CRCL 150 mL/min",
              "AST 15 U/L","AST 34 U/L",
              "ALT 14 U/L","ALT 43 U/L",
              "BILI 5 umol/L","BILI 15 umol/L",
              "Male","Female")

## The print names of the covariates
covariates <- c("NCI","Formulation","Food status","2D6 genotype","Race","Weight","BMI","Age","Createnine\nclearance","AST","ALT","Bilirubin","Sex")

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


extFile <- file.path(modDevDir ,paste0("run",finalRun,".ext")) # The path to the ext file
covFile <- file.path(modDevDir ,paste0("run",finalRun,".cov")) # The path path to the cov file
dfSamplesCOV <- getSamples(covFile,extFile=extFile,n=200)


dfres <- getForestDFSCM(dfCovs           = dfCovs,
                        cdfCovsNames     = covnames,
                        functionList     = list(paramFunction),
                        functionListName = functionListName,
                        noBaseThetas     = 14,
                        dfParameters     = dfSamplesCOV,
                        dfRefRow         = NULL
)



old <- theme_set(theme_bw(base_size = 12))

forestPlot(dfres,parameters = c("CL","Frel","AUC"),groupNameLabels=covariates,size=16,tabTextSize = 16,
           referenceInfo = "Reference subject: 85 kg, Male, FDC, Fed, CYP2D6 PM",strip.top.size = 6)


