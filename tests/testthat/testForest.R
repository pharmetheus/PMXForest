context("Forest plot functions")
library(testthat)
library(PMXForest)


test_that("dfCovs is created properly", {
  dfCovs1 <- createInputForestData(
    list( "FORM" = c(0,1),
          "FOOD" = c(0,1),
          "GENO" = c(1,2,3,4),
          "RACEL"= c(1,2,3),
          "WT"   = c(65,115),
          "AGE"  = c(27,62),
          "CRCL" = c(83,150),
          "SEX"  = c(1,2)),
    iMiss=-99)

  expect_equal_to_reference(dfCovs1,"test_output/dfCovs1Output")
})

test_that("getForestDFSCM works properly", {

  set.seed(123)
  runno   <- 1
  extFile <- system.file("extdata",paste0("run",runno,".ext"),package="PMXForest")
  covFile <- system.file("extdata",paste0("run",runno,".cov"),package="PMXForest")

  dfSamplesCOV     <- getSamples(covFile,extFile=extFile,n=200)

  paramFunction <- function(thetas, df, ...) {

    FRELNCIL <- 1 # Most common
    if (any(names(df) == "NCIL") && df$NCIL == 1) FRELNCIL <- 1 + thetas[16]

    FRELFORM <- 1 # Most common
    if (any(names(df) == "FORM") && df$FORM == 0) FRELFORM <- 1 + thetas[15]

    FRELCOV <- FRELFORM * FRELNCIL

    CLFOOD <- 1 # Most common
    if (df$FOOD == 0) CLFOOD <- 1 + thetas[14]

    CLCOV <- CLFOOD

    TVFREL <- thetas[1]
    if (df$GENO == 1) TVFREL <- TVFREL * (1 + thetas[11])
    if (df$GENO == 3) TVFREL <- TVFREL * (1 + thetas[12])
    if (df$GENO == 4) TVFREL <- TVFREL * (1 + thetas[13])

    TVFREL <- FRELCOV * TVFREL

    if (!any(names(df) == "WT") || df$WT == -99) {
      TVCL <- thetas[4]
    } else {
      TVCL <- thetas[4] * (df$WT / 75)**thetas[2]
    }

    if (df$GENO == 1) TVCL <- TVCL * (1 + thetas[8])
    if (df$GENO == 3) TVCL <- TVCL * (1 + thetas[9])
    if (df$GENO == 4) TVCL <- TVCL * (1 + thetas[10])

    TVCL <- CLCOV * TVCL

    if (!any(names(df) == "WT") || df$WT == -99) {
      TVV <- thetas[5]
    } else {
      TVV <- thetas[5] * (df$WT / 75)**thetas[3]
    }
    TVMAT <- thetas[6]
    TVD1 <- thetas[7]

    FREL <- TVFREL
    CL <- TVCL
    V <- TVV

    if (any(names(df) == "REFROW")) CL <- CL$REFROW

    return(CL)
  }

  functionListName <- c("CL")

  dfCovs1 <- createInputForestData(
    list( "FORM" = c(0,1),
          "FOOD" = c(0,1),
          "GENO" = c(1,2,3,4),
          "RACEL"= c(1,2,3),
          "WT"   = c(65,115),
          "AGE"  = c(27,62),
          "CRCL" = c(83,150),
          "SEX"  = c(1,2)),
    iMiss=-99)

  covnames <- c("Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM",
                "Caucasian","African American","Asian and other","WT 65 kg","WT 115 kg",
                "Age 27 y","Age 62 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

  dfRefRow <- data.frame("FORM"=1,"FOOD"=1,RACEL=0,"GENO"=2,"WT"=75,"AGE"=70,"CRCL"=118,"SEX"=1)

  dfresCOV <- getForestDFSCM(
    dfCovs = dfCovs1,
    cdfCovsNames = covnames,
    functionList = list(paramFunction),
    functionListName = functionListName,
    noBaseThetas = 16,
    noSigmas = 1,
    noSkipOm = 0,
    dfParameters = dfSamplesCOV,
    probs = c(0.05, 0.5, 0.95),
    dfRefRow = dfRefRow,
    quiet = TRUE,
    groupdist = 0.3,
    withingroupdist = 0.2)

  expect_equal_to_reference(dfresCOV,"test_output/dfresCovOutput")

  ## Test that the cores argument works
  dfresCOVcores <- getForestDFSCM(
    dfCovs = dfCovs1,
    cdfCovsNames = covnames,
    functionList = list(paramFunction),
    functionListName = functionListName,
    noBaseThetas = 16,
    noSigmas = 1,
    noSkipOm = 0,
    dfParameters = dfSamplesCOV,
    probs = c(0.05, 0.5, 0.95),
    dfRefRow = dfRefRow,
    quiet = TRUE,
    groupdist = 0.3,
    withingroupdist = 0.2,
    cores=2)

  expect_equal_to_reference(dfresCOVcores,"test_output/dfresCovOutput")

  p1 <- plotForestDF(dfresCOV,textsize = 5)+
    xlim(0,3.3) +
    geom_vline(xintercept=c(0.8,1.25),linetype="dashed")

  expect_equal_to_reference(p1,"test_output/plotForestDFoutput")


  ## Test that dfRefRow works
  dfRefRowTest <- data.frame("FORM"=1,"FOOD"=1,RACEL=0,"GENO"=1,"WT"=75,"AGE"=70,"CRCL"=118,"SEX"=1)

  dfresCOVdfRefRowTest1 <- getForestDFSCM(
    dfCovs = dfCovs1,
    cdfCovsNames = covnames,
    functionList = list(paramFunction),
    functionListName = functionListName,
    noBaseThetas = 16,
    noSigmas = 1,
    noSkipOm = 0,
    dfParameters = dfSamplesCOV,
    probs = c(0.05, 0.5, 0.95),
    dfRefRow = dfRefRowTest,
    quiet = TRUE,
    groupdist = 0.3,
    withingroupdist = 0.2,
    cores=2)

  expect_equal_to_reference(dfresCOVdfRefRowTest1,"test_output/dfresCovdfRefRowTest1")

  dfresCOVdfRefRowTestNULL <- getForestDFSCM(
    dfCovs = dfCovs1,
    cdfCovsNames = covnames,
    functionList = list(paramFunction),
    functionListName = functionListName,
    noBaseThetas = 16,
    noSigmas = 1,
    noSkipOm = 0,
    dfParameters = dfSamplesCOV,
    probs = c(0.05, 0.5, 0.95),
    dfRefRow = NULL,
    quiet = TRUE,
    groupdist = 0.3,
    withingroupdist = 0.2,
    cores=2)

  expect_equal_to_reference(dfresCOVdfRefRowTestNULL,"test_output/dfresCovdfRefRowTestNULL")

  ## Test the empirical forest plots

  dataFile <- system.file("extdata","DAT-1-MI-PMX-2.csv",package="PMXForest")

  dfData <- read.csv(dataFile) %>% distinct(ID,.keep_all=TRUE)
  dfData <- dfData[1:300,]

  lsExpr<-rev(list(expression(SEX==2),
                   expression(SEX==1),
                   expression(CRCL>=146),
                   expression(CRCL<94),
                   expression(AGE>=57),
                   expression(AGE<35),
                   expression(WT>104),
                   expression(WT<70),
                   expression(RACEL==3),
                   expression(RACEL==2),
                   expression(RACEL==1),
                   expression(GENO==4),
                   expression(GENO==3),
                   expression(GENO==2),
                   expression(GENO==1),
                   expression(FOOD==1),
                   expression(FOOD==0),
                   expression(FORM==1),
                   expression(FORM==0)))

  cGrouping = c(1,1,2,2,3,3,3,3,4,4,4,5,5,6,6,7,7,8,8,9,9)

  covnamesEmp <- c("Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian","African American","Asian and other","WT<70 kg","WT>104 kg",
                   "Age<35 y","Age>57 y","CRCL<94 mL/min","CRCL>146 mL/min","Male","Female")


  ## $COV uncertainty
  dfresCOVemp<-getForestDFemp(
    metricFunction     = median,
    cGrouping          = cGrouping,
    dfData             = dfData,
    covExpressionsList = lsExpr,
    cdfCovsNames       = covnamesEmp,
    functionList       = list(paramFunction),
    functionListName   = functionListName,
    noBaseThetas       = 16,
    dfParameters       = dfSamplesCOV,
    ncores             = 6)

  expect_equal_to_reference(p1,"test_output/dfresCOVempoutput")

  p2 <- plotForestDF(dfresCOVemp,textsize = 5)+
    xlim(0,3.3) +
    geom_vline(xintercept=c(0.8,1.25),linetype="dashed")

  expect_equal_to_reference(p2,"test_output/plotForestDFempoutput")

})
