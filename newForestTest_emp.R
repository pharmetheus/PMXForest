library(dplyr)
library(ggplot2)
library(ggpubr)
## This file sets up the reference objects for the ubuts tests
set.seed(865765)

lsExpr<-rev(list("SEX"=expression(SEX==2),
                 "SEX"=expression(SEX==1),
                 "CRCL"=expression(CRCL>=146),
                 "CRCL"=expression(CRCL<94),
                 "AGE"=expression(AGE>=57),
                 "AGE"=expression(AGE<35),
                 "WT"=expression(WT>104),
                 "WT"=expression(WT<70),
                 "RACE"=expression(RACEL==3),
                 "RACE"=expression(RACEL==2),
                 "RACE"=expression(RACEL==1),
                 "GENO"=expression(GENO==4),
                 "GENO"=expression(GENO==3),
                 "GENO"=expression(GENO==2),
                 "GENO"=expression(GENO==1),
                 "FOOD"=expression(FOOD==1),
                 "FOOD"=expression(FOOD==0),
                 "FORM"=expression(FORM==1),
                 "FORM"=expression(FORM==0)))

covnamesEmp <- c("Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian","African American","Asian and other","WT<70 kg","WT>104 kg",
                 "Age<35 y","Age>57 y","CRCL<94 mL/min","CRCL>146 mL/min","Male","Female")

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

dataFile <- system.file("extdata","DAT-1-MI-PMX-2.csv",package="PMXForest")
dfData   <- read.csv(dataFile) %>% distinct(ID,.keep_all=TRUE)


dfresEmp <- getForestDFemp(
  dfData             = dfData,
  covExpressionsList = lsExpr,
  cdfCovsNames       = covnamesEmp,
  functionList       = list(paramFunction),
  functionListName   = functionListName,
  metricFunction     = mean,
  noBaseThetas       = 16,
  dfParameters       = dfSamplesCOV,
  dfRefRow           = NULL,
  ncores             = 6,
  cstrPackages       = "dplyr"
)



plotData       <- setupForestPlotData(dfresEmp,plotRelative=FALSE,noVar=TRUE)
# plotDataRefRow <- setupForestPlotData(dfresRefRow,plotRelative=FALSE,noVar=TRUE,sigdig=3)


old <- theme_set(theme_bw())

forestPlot(dfresEmp,plotRelative=FALSE)
# forestPlot(dfres,tabTextSize=5)
# forestPlot(dfresRefRow,plotRelative = FALSE)
# forestPlot(dfres,plotData=plotData)
# forestPlot(dfres,plotRelative = FALSE)
# forestPlot(dfres,parameters = c("CL","FREL"),parameterLabels=c("CL(L/h)","F[rel]"),groupNameLabels=covariates)
# forestPlot(dfres,plotData=myPlotData,parameters = c("CL","FREL"),parameterLabels=c("CL(L/h)","F[rel]"),groupNameLabels=covariates)
