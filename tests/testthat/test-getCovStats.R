test_that("getCovStats works", {

  library(PMXFrem)
  modDevDir <- system.file("extdata/SimNeb",package = "PMXFrem")
  fremRun   <- 31
  fremFiles <- getFileNames(31,modDevDir = modDevDir)
  covNames  <- getCovNames(fremFiles$mod)

  data <- fread(file.path(modDevDir,"DAT-2-MI-PMX-2-onlyTYPE2-new.csv"))

  dfCovs <- PMXForest::createInputForestData(
    getCovStats(data,covNames$orgCovNames,missVal=-99),
    iMiss    =-99
  )

  expect_snapshot(dfCovs)

  dfCovs2 <- PMXForest::createInputForestData(
    getCovStats(data,covNames$orgCovNames,missVal=-99,probs=c(0.1,0.9)),
    iMiss    =-99
  )

  expect_snapshot(dfCovs2)
})
