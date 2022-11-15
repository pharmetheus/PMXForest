library(testthat)
library("PMXFrem",lib.loc="~/4.1/library/")

test_that("Forest plots for FREM works properly", {

#  if(any(.packages(all.available = T,lib.loc="~/4.1/library/")=="PMXFrem")) skip("PMXFrem not available")

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


  expect_equal_to_reference(dfCovs,"test_output/dfCovsOutputFREM")

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

  set.seed(123)
  runno   <- "22-3"

  extFile <- system.file("extdata",paste0("SimVal/run",runno,".ext"),package="PMXForest")
  covFile <- system.file("extdata",paste0("SimVal/run",runno,".cov"),package="PMXForest")
  modFile <- system.file("extdata",paste0("SimVal/run",runno,".mod"),package="PMXForest")

  ## Will use only 25 samples for the tests
  dfSamplesCOV     <- getSamples(covFile,extFile=extFile,n=25)

  covnames  <- c("NCI=0","NCI>0","Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian",
                 "Other","WT 64 kg","WT 116 kg",
                 "BMI 23 kg/m^2","BMI 39 kg/m^2",
                 "Age 25 y","Age 61 y","CRCL 83 mL/min","CRCL 150 mL/min","AST 15 IU","AST 34 IU","ALT 14 IU","ALT 43 IU",
                 "Bilirubin 5 µmol/L", "Bilirubin 15 µmol/L","Male","Female")

  covariates <- c("NCI","Formulation","Food status","2D6 genotype","Race","Weight","BMI","Age","Createnine\nclearance","AST","ALT","Bilirubin","Sex")

  numNonFREMThetas <- 13
  numSkipOm <- 2

  dfresFREM <-getForestDFFREM(dfCovs = dfCovs,
                              cdfCovsNames = covnames,
                              covNames=PMXFrem::getCovNames(modFile),
                              functionList = list(paramFunction),
                              functionListName = functionListName,
                              numSkipOm =numSkipOm,
                              numNonFREMThetas =numNonFREMThetas,
                              # noBaseThetas = noBaseThetas,
                              # noParCov = 4,
                              # noSigmas = 2,
                              # noSkipOm = 2,
                              # noCovThetas = noCovThetas,
                              dfParameters = dfSamplesCOV,
                              probs = c(0.05, 0.95),
                              dfRefRow = NULL,
                              quiet = TRUE,
                              ncores = 1,
                              cstrPackages = c("PMXFrem","dplyr"))

  dfresFREM2 <-getForestDFFREM(dfCovs          = dfCovs,
                               cdfCovsNames     = covnames,
                               covNames         = PMXFrem::getCovNames(modFile),
                               functionList     = list(paramFunction),
                               functionListName = functionListName,
                               numSkipOm =numSkipOm,
                               numNonFREMThetas =numNonFREMThetas,
                               # noBaseThetas     = noBaseThetas,
                               # noParCov         = 4,
                               # noSigmas         = 2,
                               # noSkipOm         = 2,
                               # noCovThetas      = noCovThetas,
                               dfParameters     = dfSamplesCOV,
                               probs            = c(0.05, 0.95),
                               dfRefRow         = dfCovs %>% filter(COVARIATEGROUPS == "BMI") %>% select(-COVARIATEGROUPS) %>% slice(1),
                               quiet            = TRUE,
                               ncores           = 1,
                               cstrPackages     = c("PMXFrem","dplyr"))


  dfresFREM3 <-getForestDFFREM(dfCovs          = dfCovs,
                               #cdfCovsNames     = covnames,
                               covNames         = PMXFrem::getCovNames(modFile),
                               functionList     = list(paramFunction),
                               functionListName = functionListName,
                               numSkipOm =numSkipOm,
                               numNonFREMThetas =numNonFREMThetas,
                               # noBaseThetas     = noBaseThetas,
                               # noParCov         = 4,
                               # noSigmas         = 2,
                               # noSkipOm         = 2,
                               # noCovThetas      = noCovThetas,
                               dfParameters     = dfSamplesCOV,
                               probs            = c(0.05, 0.95),
                               quiet            = TRUE,
                               ncores           = 1,
                               cstrPackages     = c("PMXFrem","dplyr"))

  # Covariates as a list instead of data.frame
  covList <- list("NCIL" = c(0,1),
                  "FORM" = c(0,1),
                  "FOOD" = c(0,1),
                  "GENO"=list("GENO1" = c(1,0,0,0),
                              "GENO3" = c(0,0,1,0),
                              "GENO4" = c(0,0,0,1)),
                  "RACEL"= c(1,2,3),
                  "WT"   = c(65,115),
                  "AGE"  = c(27,62),
                  "CRCL" = c(83,150),
                  "SEX"  = c(1,2))

  dfresFREM4 <-getForestDFFREM(dfCovs           = covList,
                               #cdfCovsNames   = covnames,
                               covNames         = PMXFrem::getCovNames(modFile),
                               functionList     = list(paramFunction),
                               functionListName = functionListName,
                               numSkipOm =numSkipOm,
                               numNonFREMThetas =numNonFREMThetas,
                               # noBaseThetas     = noBaseThetas,
                               # noParCov         = 4,
                               # noSigmas         = 2,
                               # noSkipOm         = 2,
                               # noCovThetas      = noCovThetas,
                               dfParameters     = dfSamplesCOV,
                               probs            = c(0.05, 0.95),
                               quiet            = TRUE,
                               ncores           = 1,
                               cstrPackages     = c("PMXFrem","dplyr"))

  rVersion <- paste0(R.version$major,".",R.version$minor)

  expect_equal_to_reference(dfresFREM, paste0("test_output/dfresFREM",rVersion))
  expect_equal_to_reference(dfresFREM2,paste0("test_output/dfresFREM2",rVersion))
  expect_equal_to_reference(dfresFREM3,paste0("test_output/dfresFREM3",rVersion))
  expect_equal_to_reference(dfresFREM4,paste0("test_output/dfresFREM4",rVersion))


  ## Tests for setupData
  plotData1 <- setupForestPlotData(dfresFREM)
  plotData2 <- setupForestPlotData(dfresFREM,plotRelative = FALSE,noVar=TRUE)
  plotData4 <- setupForestPlotData(dfresFREM,parameterLabels=c("CL (L/h)","V (L)","Frel"),groupNameLabels=covariates,plotRelative=FALSE,noVar=TRUE,statisticsLabel=c("CL (L/h)","V (L)","Frel"))
  expect_equal_to_reference(plotData1,"test_output/plotData1frem")
  expect_equal_to_reference(plotData2,"test_output/plotData2frem")
  expect_equal_to_reference(plotData4,"test_output/plotData4frem")

  check_graphical_output <- function() {
    if(getRversion() <= "3.5.3") {
      skip("R version <= 3.5.3")
    } else {

      svg() # Start a device to make plots cosistent between different ways of running the tests
      fp5 <- forestPlot(dfresFREM)
      fp6 <- forestPlot(dfresFREM,plotData=plotData1)
      fp7 <- forestPlot(dfresFREM,plotData=plotData1 %>%
                          filter(PARAMETER=="CL") %>%
                          group_by(GROUPNAME) %>%
                          mutate(COVEFF=ifelse(any(COVEFF==TRUE),TRUE,FALSE)) %>%
                          filter(COVEFF==TRUE),parameters="CL")
      fp8 <- forestPlot(dfresFREM,plotRelative = FALSE)
      fp9 <- forestPlot(dfresFREM,parameters="CL")
      fp10 <- forestPlot(dfresFREM,rightStrip = FALSE)
      fp11 <- forestPlot(dfresFREM,rightStrip = FALSE,table=FALSE)
      fp12 <- forestPlot(dfresFREM,table=FALSE)

      dev.off()

      vdiffr::expect_doppelganger("Forest plot with default options", fp5)
      vdiffr::expect_doppelganger("Forest plot with provided plot data", fp6)
      vdiffr::expect_doppelganger("Forest plot with subsetting", fp7)
      vdiffr::expect_doppelganger("Forest plot with plotRelative=FALS", fp8)
      vdiffr::expect_doppelganger("Forest plot only CL", fp9)
      vdiffr::expect_doppelganger("Forest plot without rightStrip", fp10)
      vdiffr::expect_doppelganger("Forest plot without rightStrip and no table", fp11)
      vdiffr::expect_doppelganger("Forest plot without table", fp12)
    }
  }

  check_graphical_output()

})


