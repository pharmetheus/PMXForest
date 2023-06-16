## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
packageName <- "PMXForest"
suppressPackageStartupMessages(library(packageName,character.only = TRUE))
library(tools)
library(ggplot2)
library(stats)



theme_set(theme_bw(base_size=18))
theme_update(plot.title = element_text(hjust = 0.5))

set.seed(865765)

## -----------------------------------------------------------------------------
dfCovs <- createInputForestData(
  list( "AGE" = c(35,50,65),
        "EXP" = c(0,100,200,300,400,500)),
  iMiss=-99) 

## ---- echo=F------------------------------------------------------------------
dfCovs

## -----------------------------------------------------------------------------
covnames <- c("Age 35 y","Age 50 y","Age 65 y","Placebo","Drug conc of 100 ng/mL","Drug conc of 200 ng/mL","Drug conc of 300 ng/mL","Drug conc of 400 ng/mL","Drug conc of 500 ng/mL")

## -----------------------------------------------------------------------------
dfRefRow <- data.frame("AGE"=50,"EXP"=0)

## -----------------------------------------------------------------------------

#Hazard ratio (compare to the reference subject) parameter function
paramFunctionHR <- function(thetas, df, ...) {

  # Risk factors / covariates section
  RF<-0 
  
  if(any(names(df) == "AGE") && df$AGE!=-99) {
    RF<-RF+(df$AGE-50)*thetas[3]
  }

  # Drug /exposure effect section
  
  EFF<-0
  if(any(names(df) == "EXP") && df$EXP!=-99) {
    EFF<-thetas[4]*df$EXP
  }

  #Hazard ratio given the covariates and exposure
  HR<-exp(RF+EFF)
  
  return(HR)
}

#Probability of having an event up until time=TIME
paramFunctionCDF <- function(thetas, df, functHR,TIME, ...) {

  #The HR (exponentiated covariate and exposure effects)
  COV<-functHR(thetas,df,...)
  
  #Re-parameterize the model to fit the pweibull parameterization
  #Note that this is dependent on the parameterization used in the model
  TH1     <- 1/thetas[1]
  TH2     <- thetas[2]

  LAM     <- ((1/COV)^(1/TH2))*TH1

  
  #Get the probability of having an event up until time TIME (CDF)
  F       <- pweibull(TIME,TH2,LAM)
  return (F)
}



## -----------------------------------------------------------------------------
functionListNameHR <- c("Hazard Ratio")
functionListNameCDF <- c("Probability of an event within 1 year")

## -----------------------------------------------------------------------------
runno   <- "tte_weibull"

extFile <- system.file("extdata","tte",paste0(runno,".ext"),package=packageName)

## -----------------------------------------------------------------------------
covFile <- system.file("extdata","tte",paste0(runno,".cov"),package=packageName) 
dfSamplesCOV <- getSamples(covFile,extFile=extFile,n=200)

## -----------------------------------------------------------------------------
bootFile    <- system.file("extdata","tte","bootstrap_tte_weibull_n500",paste0("raw_results_",runno,".csv"),package=packageName)
dfSamplesBS <- getSamples(bootFile,extFile=extFile)

## -----------------------------------------------------------------------------
#Hazard ratio forest data
  dfresHR <- getForestDFSCM(dfCovs = dfCovs,
                          cdfCovsNames = covnames,
                          functionList = list(paramFunctionHR),
                          functionListName = functionListNameHR,
                          noBaseThetas = 4,
                          dfParameters = dfSamplesBS, 
                          dfRefRow = dfRefRow,
                          cstrPackages = c("dplyr"))

#Probability of event within one year forest data
  dfresCDF <- getForestDFSCM(dfCovs = dfCovs,
                          cdfCovsNames = covnames,
                          functionList = list(paramFunctionCDF),
                          functionListName = functionListNameCDF,
                          noBaseThetas = 4,
                          dfParameters = dfSamplesBS, 
                          dfRefRow = dfRefRow,
                          cstrPackages = c("dplyr","stats"),
                          functHR=paramFunctionHR,
                          TIME=365)



## ----ForestHR,fig.width=12,fig.asp=.5,out.width="100%"------------------------

plotList<-forestPlot(dfresHR,plotRelative = FALSE,xlb = "Hazard Ratio",return="plotList")
#Add the log10 x-axis
plotList[[1]]<-plotList[[1]]+scale_x_continuous(trans='log10')

ggpubr::ggarrange(plotlist      = plotList,
                  ncol          = 2,
                  nrow          = 1,
                  widths        = c(2,0.7),
                  align         ="h",
                  common.legend = T)

## ----ForestCDF,fig.width=12,fig.asp=0.6,out.width="100%"----------------------


plotList<-forestPlot(dfresCDF,plotRelative = FALSE,xlb = "Probability of event",return="plotList")
#Add the Percent x-axis
plotList[[1]]<-plotList[[1]]+scale_x_continuous(labels = scales::percent)


ggpubr::ggarrange(plotlist      = plotList,
                  ncol          = 2,
                  nrow          = 1,
                  widths        = c(1,0.7),
                  align         ="h",
                  common.legend = T)


