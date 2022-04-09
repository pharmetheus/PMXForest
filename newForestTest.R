library(dplyr)
library(ggplot2)
library(ggpubr)
## This file sets up the reference objects for the ubuts tests
set.seed(865765)

dfCovs <- createInputForestData(
  list( "FORM" = c(0,1),
        "FOOD" = c(0,1),
        "GENO" = c(1,2,3,4),
        "RACEL"= c(1,2,3),
        "WT"   = c(65,115),
        "AGE"  = c(27,62),
        "CRCL" = c(83,150),
        "SEX"  = c(1,2)),
  iMiss=-99)

#cGrouping <- c(1,1,2,2,3,3,3,3,4,4,4,5,5,6,6,7,7,8,8)
covnames  <- c("Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian",
               "African American","Asian and other","WT 65 kg","WT 115 kg",
               "Age 27 y","Age 62 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

## The print names of the covariates
covariates <- c("Formulation","Food status","2D6 genotype","Race","Weight","Age","Createnine\nclearance","Sex")

paramFunction <- function(thetas, df, ...) {

  FRELNCIL <- 1
  if (any(names(df) == "NCIL") && df$NCIL == 1) FRELNCIL <- 1 + thetas[16]

  FRELFORM <- 1
  if (any(names(df) == "FORM") && df$FORM == 0) FRELFORM <- 1 + thetas[15]

  FRELCOV <- FRELFORM * FRELNCIL

  CLFOOD <- 1
  if (df$FOOD == 0) CLFOOD <- 1 + thetas[14]

  CLCOV <- CLFOOD

  TVFREL <- thetas[1]
  if(any(names(df) == "GENO")) {
    if(df$GENO == 1) TVFREL <- TVFREL * (1 + thetas[11])
    if(df$GENO == 3) TVFREL <- TVFREL * (1 + thetas[12])
    if(df$GENO == 4) TVFREL <- TVFREL * (1 + thetas[13])
  }

  TVFREL <- FRELCOV * TVFREL

  if (!any(names(df) == "WT") || df$WT == -99) {
    TVCL <- thetas[4]
  } else {
    TVCL <- thetas[4] * (df$WT / 75)**thetas[2]
  }

  if(any(names(df) == "GENO")) {
    if (df$GENO == 1) TVCL <- TVCL * (1 + thetas[8])
    if (df$GENO == 3) TVCL <- TVCL * (1 + thetas[9])
    if (df$GENO == 4) TVCL <- TVCL * (1 + thetas[10])
  }

  TVCL <- CLCOV * TVCL

  if (!any(names(df) == "WT") || df$WT == -99) {
    TVV <- thetas[5]
  } else {
    TVV <- thetas[5] * (df$WT / 75)**thetas[3]
  }
  TVMAT <- thetas[6]
  TVD1  <- thetas[7]

  FREL <- TVFREL
  CL <- TVCL
  V <- TVV

  return(list(CL,FREL,V))
}

functionListName <- c("CL","FREL","V")


runno   <- 1
extFile <- system.file("extdata",paste0("run",runno,".ext"),package="PMXForest")
covFile <- system.file("extdata",paste0("run",runno,".cov"),package="PMXForest")
dfSamplesCOV <- getSamples(covFile,extFile=extFile,n=200)


dfres <- getForestDFSCM(dfCovs           = dfCovs,
                        cdfCovsNames     = covnames,
                        functionList     = list(paramFunction),
                        functionListName = functionListName,
                        noBaseThetas     = noBaseThetas,
                        dfParameters     = dfSamplesCOV,
                        dfRefRow         = NULL
)


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
forestPlot(dfresRefRow,plotRelative = FALSE)
forestPlot(dfres,plotData=plotData)
forestPlot(dfres,plotRelative = FALSE)
forestPlot(dfres,parameters = c("CL","FREL"),parameterLabels=c("CL(L/h)","F[rel]"),groupNameLabels=covariates)
forestPlot(dfres,plotData=myPlotData,parameters = c("CL","FREL"),parameterLabels=c("CL(L/h)","F[rel]"),groupNameLabels=covariates)
