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
  dfSamplesCOV     <- getSamples(covFile,extFile=extFile,n=175)

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


  expect_equal_to_reference(dfresSCM,"test_output/dfresSCM")


  ## Tests for setupData
  plotData1 <- setupForestPlotData(dfresSCM)
  plotData2 <- setupForestPlotData(dfresSCM,plotRelative = FALSE,noVar=TRUE)
  plotData4 <- setupForestPlotData(dfresSCM,parameterLabels=c("CL (L/h)","V (L)","Frel"),groupNameLabels=covariates,plotRelative=FALSE,noVar=TRUE,statisticsLabel=c("CL (L/h)","V (L)","Frel"))
  expect_equal_to_reference(plotData1,"test_output/plotData1")
  expect_equal_to_reference(plotData2,"test_output/plotData2")
  expect_equal_to_reference(plotData4,"test_output/plotData4")


  ## Create a number of Forest plots to test functionality
  # fp5 <- forestPlot(dfresSCM)
  # fp6 <- forestPlot(dfresSCM,plotData=plotData1)
  # fp7 <- forestPlot(dfresSCM,plotData=plotData1 %>%
  #                     filter(PARAMETER=="CL") %>%
  #                     group_by(GROUPNAME) %>%
  #                     mutate(COVEFF=ifelse(any(COVEFF==TRUE),TRUE,FALSE)) %>%
  #                     filter(COVEFF==TRUE),parameters="CL")
  # fp8 <- forestPlot(dfresSCM,plotRelative = FALSE)
  # fp9 <- forestPlot(dfresSCM,parameters="CL")
  # fp10 <- forestPlot(dfresSCM,rightStrip = FALSE)
  # fp11 <- forestPlot(dfresSCM,rightStrip = FALSE,table=FALSE)
  # fp12 <- forestPlot(dfresSCM,table=FALSE)
  # fp13 <- forestPlot(dfresSCM,parameters = c("CL","Frel"),
  #                    stackedPlots = TRUE,
  #                    keepYlabs = TRUE,
  #                    keepRightStrip = TRUE)

  # expect_known_value(fp5,"test_output/fp5")
  # expect_equal_to_reference(fp6,"test_output/fp6")
  # expect_equal_to_reference(fp7,"test_output/fp7")
  # expect_equal_to_reference(fp8,"test_output/fp8")
  # expect_equal_to_reference(fp9,"test_output/fp9")
  # expect_equal_to_reference(fp10,"test_output/fp10")
  # expect_equal_to_reference(fp11,"test_output/fp11")
  # expect_equal_to_reference(fp12,"test_output/fp12")
  # expect_equal_to_reference(fp13,"test_output/fp13")
})


