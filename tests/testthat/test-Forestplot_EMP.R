dontCheckGraphs <- FALSE

test_that("Forest plots for EMP works properly", {

  # skip("Save run time")

  lsExpr<-rev(
    list("SEX" = expression(SEX==2),
         "SEX" = expression(SEX==1),
         "CRCL"= expression(CRCL>=146),
         "CRCL"= expression(CRCL<94),
         "AGE" = expression(AGE>=57),
         "AGE" = expression(AGE<35),
         "WT"  = expression(WT>104),
         "WT"  = expression(WT<70),
         "RACE"= expression(RACEL==3),
         "RACE"= expression(RACEL==2),
         "RACE"= expression(RACEL==1),
         "GENO"= expression(GENO==4),
         "GENO"= expression(GENO==3),
         "GENO"= expression(GENO==2),
         "GENO"= expression(GENO==1),
         "FOOD"= expression(FOOD==1),
         "FOOD"= expression(FOOD==0),
         "FORM"= expression(FORM==1),
         "FORM"= expression(FORM==0),
         "NCI" = expression(NCIL==1),
         "NCI" = expression(NCIL==0)
         )
  )


  expect_snapshot(lsExpr)

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
  functionList<-list(paramFunction)

  set.seed(123)
  runno   <- 7
  extFile <- system.file("extdata",paste0("SimVal/run",runno,".ext"),package="PMXForest")
  covFile <- system.file("extdata",paste0("SimVal/run",runno,".cov"),package="PMXForest")

  ## Will use only 25 samples for the tests
  dfSamplesCOV     <- getSamples(covFile,extFile=extFile,n=25)

  dataFile <- system.file("extdata","SimVal/DAT-1-MI-PMX-2.csv",package="PMXForest")
  dfData   <- read.csv(dataFile,stringsAsFactors = FALSE) %>% distinct(ID,.keep_all=TRUE) %>% slice(1:20,500:600)

  noBaseThetas <- 14

  covnamesEmp <- c("NCI=0","NCI>0","Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian","African American","Asian and other","WT<70 kg","WT>104 kg",
                   "Age<35 y","Age>57 y","CRCL<94 mL/min","CRCL>146 mL/min","Male","Female")

  covariates <- c("Formulation","Food status","2D6 genotype","Race","Weight","Age","Createnine\nclearance","Sex","NCI")

  dfresEMP <- getForestDFemp(dfData             = dfData,
                             covExpressionsList = lsExpr,
                             cdfCovsNames       = covnamesEmp,
                             functionList       = list(paramFunction),
                             functionListName   = functionListName,
                             metricFunction     = median,
                             noBaseThetas       = noBaseThetas,
                             dfParameters       = dfSamplesCOV,
                             dfRefRow           = NULL,
                             ncores             = 1,
                             cstrPackages       = "dplyr"
  )

  lsExprRef<-rev(
    list("SEX" = expression(SEX==2))
    )

  dfresEMP2 <- getForestDFemp(dfData             = dfData,
                              covExpressionsList = lsExpr,
                              cdfCovsNames       = covnamesEmp,
                              functionList       = list(paramFunction),
                              functionListName   = functionListName,
                              metricFunction     = median,
                              noBaseThetas       = noBaseThetas,
                              dfParameters       = dfSamplesCOV,
                              dfRefRow           = lsExprRef,
                              ncores             = 1,
                              cstrPackages       = "dplyr"
  )

  rVersion <- paste0(R.version$major,".",R.version$minor)
  expect_snapshot(dfresEMP  )
  expect_snapshot(dfresEMP2 )


  ## Tests for setupData
  plotData1 <- setupForestPlotData(dfresEMP)
  plotData2 <- setupForestPlotData(dfresEMP,plotRelative = FALSE,noVar=TRUE)
  plotData4 <- setupForestPlotData(dfresEMP,parameterLabels=c("CL (L/h)","V (L)","Frel"),groupNameLabels=covariates,plotRelative=FALSE,noVar=TRUE,statisticsLabel=c("CL (L/h)","V (L)","Frel"))
  expect_snapshot(plotData1)
  expect_snapshot(plotData2)
  expect_snapshot(plotData4)

  check_graphical_output <- function() {
    if(getRversion() <= "3.5.3") {
      skip("R version <= 3.5.3")
    } else {

      if(dontCheckGraphs) {
        skip("Don't check plots")
      }

      svg() # Start a device to make plots cosistent between different ways of running the tests
      fp5 <- forestPlot(dfresEMP)
      fp6 <- forestPlot(dfresEMP,plotData=plotData1)
      fp7 <- forestPlot(dfresEMP,plotData=plotData1 %>%
                          filter(PARAMETER=="CL") %>%
                          group_by(GROUPNAME) %>%
                          mutate(COVEFF=ifelse(any(COVEFF==TRUE),TRUE,FALSE)) %>%
                          filter(COVEFF==TRUE),parameters="CL")
      fp8 <- forestPlot(dfresEMP,plotRelative = FALSE)
      fp9 <- forestPlot(dfresEMP,parameters="CL")
      fp10 <- forestPlot(dfresEMP,rightStrip = FALSE)
      fp11 <- forestPlot(dfresEMP,rightStrip = FALSE,table=FALSE)
      fp12 <- forestPlot(dfresEMP,table=FALSE)

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


