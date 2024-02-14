test_that("getForestDFFREM works", {

  runno   <- "22-3"

  extFile <- system.file("extdata",paste0("SimVal/run",runno,".ext"),package="PMXForest")
  covFile <- system.file("extdata",paste0("SimVal/run",runno,".cov"),package="PMXForest")
  modFile <- system.file("extdata",paste0("SimVal/run",runno,".mod"),package="PMXForest")
  datFile <- system.file("extdata",paste0("SimVal/DAT-1-MI-PMX-2.csv"),package="PMXForest")

  dfData <- read.csv(datFile)

  covNames  <- PMXFrem::getCovNames(modFile = modFile)

  dfCovs <- createInputForestData(getCovStats(dfData,covNames$orgCovNames,probs=c(0.05,0.95)))

  paramFun <- function(basethetas,covthetas, dfrow, ...) {
    CL <- basethetas[1]*exp(covthetas[1])
    V <-  basethetas[2]*exp(covthetas[2])
    AUC <- 5/CL

    return(c(CL, V, AUC))
  }

  functionListName2 <- c("CL","V","AUC")

  set.seed(123)
  ## Will use only 25 samples for the tests
  dfSamplesCOV     <- getSamples(covFile,extFile=extFile,n=25)

  dfresFREM <-getForestDFFREM(dfCovs           = dfCovs,
                              covNames         = covNames$covNames,
                              functionList     = list(paramFun),
                              functionListName = functionListName2,
                              numNonFREMThetas = 13,
                              numSkipOm        = 2,
                              dfParameters     = dfSamplesCOV,
                              probs            = c(0.05, 0.95),
                              dfRefRow         = NULL,
                              quiet            = TRUE,
                              ncores           = 1,
                              cstrPackages     = c("PMXFrem","dplyr"))

  expect_snapshot(dfresFREM)

  covlabels  <- c("Age 25 y","Age 61 y","ALT 14 IU","ALT 43 IU", "AST 15 IU","AST 34 IU",
                 "Bilirubin 5 µmol/L", "Bilirubin 15 µmol/L", "BMI 23 kg/m^2","BMI 39 kg/m^2",
                 "CRCL 83 mL/min","CRCL 150 mL/min", "Caucasian","Other",
                 "HT 152 cm","HT 185 cm",
                 "NCI=0","NCI>0","White","Other",
                 "Male","Female")

  dfresFREM2 <-getForestDFFREM(dfCovs           = dfCovs,
                               cdfCovsNames     = covlabels,
                               covNames         = covNames$covNames,
                               functionList     = list(paramFun),
                               functionListName = functionListName2,
                               numNonFREMThetas = 13,
                               numSkipOm        = 2,
                               dfParameters     = dfSamplesCOV,
                               probs            = c(0.05, 0.95),
                               dfRefRow         = NULL,
                               quiet            = TRUE,
                               ncores           = 1,
                               cstrPackages     = c("PMXFrem","dplyr"))

  expect_snapshot(dfresFREM2)
})
