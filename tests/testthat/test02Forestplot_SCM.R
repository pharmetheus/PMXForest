library(testthat)

test_that("Forest plots for SCM works properly", {



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


  expect_equal_to_reference(dfCovs,"test_output/dfCovsOutput")

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

    return(list(CL,FREL,AUC))
  }

  functionListName <- c("CL","Frel","AUC")

  set.seed(123)
  runno   <- 7
  extFile <- system.file("extdata",paste0("SimVal/run",runno,".ext"),package="PMXForest")
  covFile <- system.file("extdata",paste0("SimVal/run",runno,".cov"),package="PMXForest")

  ## Will use only 25 samples for the tests
  dfSamplesCOV     <- getSamples(covFile,extFile=extFile,n=25)

  covnames  <- c("NCI=0","NCI>0","Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian",
                 "African American","Asian and other","WT 65 kg","WT 115 kg",
                 "Age 27 y","Age 62 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

  ## The print names of the covariates
  covariates <- c("NCI","Formulation","Food status","2D6 genotype","Race","Weight","Age","Createnine\nclearance","Sex")

  noBaseThetas <- 14

  dfresSCM <- getForestDFSCM(dfCovs           = dfCovs,
                             cdfCovsNames     = covnames,
                             functionList     = list(paramFunction),
                             functionListName = functionListName,
                             noBaseThetas     = noBaseThetas,
                             dfParameters     = dfSamplesCOV,
                             dfRefRow         = NULL
  )

  ## Check some variants of getForestDFSCM

  # RefRow
  dfresSCM2 <- getForestDFSCM(dfCovs           = dfCovs,
                             cdfCovsNames     = covnames,
                             functionList     = list(paramFunction),
                             functionListName = functionListName,
                             noBaseThetas     = noBaseThetas,
                             dfParameters     = dfSamplesCOV,
                             dfRefRow         = dfCovs %>% filter(COVARIATEGROUPS == "GENO") %>% slice(1)
  )

  # No covnames
  dfresSCM3 <- getForestDFSCM(dfCovs          = dfCovs,
                             # cdfCovsNames     = covnames,
                             functionList     = list(paramFunction),
                             functionListName = functionListName,
                             noBaseThetas     = noBaseThetas,
                             dfParameters     = dfSamplesCOV,
                             dfRefRow         = NULL
  )

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

  dfresSCM4 <- getForestDFSCM(dfCovs          = covList,
                              # cdfCovsNames     = covnames,
                              functionList     = list(paramFunction),
                              functionListName = functionListName,
                              noBaseThetas     = noBaseThetas,
                              dfParameters     = dfSamplesCOV,
                              dfRefRow         = NULL
  )

  rVersion <- paste0(R.version$major,".",R.version$minor)
  expect_equal_to_reference(dfresSCM,paste0("test_output/dfresSCM",rVersion))
  expect_equal_to_reference(dfresSCM3,paste0("test_output/dfresSCM3",rVersion))
  expect_equal_to_reference(dfresSCM4,paste0("test_output/dfresSCM4",rVersion))

  ## Tests for setupData
  expect_error(setupForestPlotData(dfresSCM,parameterLabels=c("CL (L/h)","V (L)","Frel","Fake")))
  expect_error(setupForestPlotData(dfresSCM,groupNameLabels=covariates[-1]))
  expect_error(setupForestPlotData(dfresSCM,statisticsLabels=c("CL (L/h)","V (L)","Frel","Fake")))


  plotData1 <- setupForestPlotData(dfresSCM)
  plotData2 <- setupForestPlotData(dfresSCM,plotRelative = FALSE,noVar=TRUE)
  plotData4 <- setupForestPlotData(dfresSCM,parameterLabels=c("CL (L/h)","V (L)","Frel"),
                                   groupNameLabels=covariates,plotRelative=FALSE,noVar=TRUE,
                                   statisticsLabel=c("CL (L/h)","V (L)","Frel"))
  plotData5 <- setupForestPlotData(dfresSCM,onlySignificant = TRUE)

  expect_equal_to_reference(plotData1,"test_output/plotData1")
  expect_equal_to_reference(plotData2,"test_output/plotData2")
  expect_equal_to_reference(plotData4,"test_output/plotData4")
  expect_equal_to_reference(plotData5,"test_output/plotData5")


  check_graphical_output <- function() {
    if(getRversion() <= "3.5.3") {
      skip("R version <= 3.5.3")
    } else {

      svg() # Start a device to make plots cosistent between different ways of running the tests
      fp5 <- forestPlot(dfresSCM)
      fp6 <- forestPlot(dfresSCM,plotData=plotData1)
      fp7 <- forestPlot(dfresSCM,plotData=plotData1 %>%
                          filter(PARAMETER=="CL") %>%
                          group_by(GROUPNAME) %>%
                          mutate(COVEFF=ifelse(any(COVEFF==TRUE),TRUE,FALSE)) %>%
                          filter(COVEFF==TRUE),parameters="CL")
      fp8 <- forestPlot(dfresSCM,plotRelative = FALSE)
      fp9 <- forestPlot(dfresSCM,parameters="CL")
      fp10 <- forestPlot(dfresSCM,rightStrip = FALSE)
      fp11 <- forestPlot(dfresSCM,rightStrip = FALSE,table=FALSE)
      fp12 <- forestPlot(dfresSCM,table=FALSE)
      fp13 <- forestPlot(dfresSCM,parameters = c("CL","Frel"),
                         stackedPlots = TRUE,
                         keepYlabs = TRUE,
                         keepRightStrip = TRUE)
      fp14 <- forestPlot(dfresSCM,parameters = c("CL","Frel"),
                         stackedPlots = TRUE,
                         keepYlabs = TRUE,
                         table     = FALSE,
                         keepRightStrip = TRUE)

      fp15 <- forestPlot(dfresSCM,referenceInfo=NULL)
      fp16 <- forestPlot(dfresSCM,referenceInfo = "Tesing of reference info)")
      fp17 <- forestPlot(dfresSCM,referenceInfo=NULL,referenceParameters = "final",noVar=FALSE)

      fp18 <- forestPlot(dfresSCM2)

      dev.off()
      vdiffr::expect_doppelganger("Forest plot with default options", fp5)
      vdiffr::expect_doppelganger("Forest plot with provided plot data", fp6)
      vdiffr::expect_doppelganger("Forest plot with subsetting", fp7)
      vdiffr::expect_doppelganger("Forest plot with plotRelative=FALS", fp8)
      vdiffr::expect_doppelganger("Forest plot only CL", fp9)
      vdiffr::expect_doppelganger("Forest plot without rightStrip", fp10)
      vdiffr::expect_doppelganger("Forest plot without rightStrip and no table", fp11)
      vdiffr::expect_doppelganger("Forest plot without table", fp12)
      vdiffr::expect_doppelganger("Forest plot stackedPlots", fp13)
      vdiffr::expect_doppelganger("Forest plot stackedPlots without table", fp14)

      vdiffr::expect_doppelganger("Forest plot referenceInfo NULL", fp15)
      vdiffr::expect_doppelganger("Forest plot referenceInfo not NULL", fp16)
      vdiffr::expect_doppelganger("Forest plot referenceParameters=final", fp17)
      vdiffr::expect_doppelganger("Forest plot stackedPlots with refRow", fp18)

    }
  }

  check_graphical_output()

})


