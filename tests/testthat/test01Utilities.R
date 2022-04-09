library(testthat)


test_that("Ext files are read properly", {
  extData <- getExt(system.file("extdata","run1.ext",package="PMXForest"))
  expect_is(extData,"data.frame")
})


test_that("mvrnorm_vector samples correctly", {

  ## Read the ext-file
  extFile <- getExt(system.file("extdata","run1.ext",package="PMXForest"))

  ## Read the cov-file and extract the covmatrix
  dfcov     <- read.table(system.file("extdata","run1.cov",package="PMXForest"),
                          fill = TRUE, header = TRUE, sep = "",
                          skip = 1, stringsAsFactors = FALSE)
  sigma     <- data.matrix(dfcov[, 2:ncol(dfcov)])

  ## Get the final parameter estimates
  finPar  <- subset(extFile, ITERATION == "-1000000000")
  mu      <- as.numeric(finPar[, -(c(1, ncol(finPar)))])

  #Get the parameters which are fixed based on the covariance
  fixedmu = rep(FALSE,1,ncol(sigma))
  for (j in 1:ncol(sigma)) {
    fixedmu[j]<-all(sigma[,j]==0)
  }

  ## Test iSampleIndex
  dfParameters <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = 10,fixed_mu=fixedmu))
  expect_equal(nrow(dfParameters),10)

  dfParameters <- mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = 1,fixed_mu = fixedmu)
  expect_is(dfParameters,"numeric")

})

test_that("getSamples works correctly for .cov input", {

  ## .cov file input
  covFile <- system.file("extdata","run1.cov",package="PMXForest")
  extFile <- system.file("extdata","run1.ext",package="PMXForest")

  expect_error(getSamples(10),"input needs to be a character string")
  expect_error(getSamples("fileNameWithoutExtension"),"The input file name needs to have an extension")
  expect_error(getSamples("fileNameWithOtherExtension.txt"),"The input file name needs to have either a .cov or .csv extension")
  expect_error(getSamples("testFile.cov"),"Can not find testFile.cov")

  expect_error(getSamples(covFile,extFile="test","test does not have a .ext extension"))
  expect_error(getSamples(covFile,"testFile.ext"),"Can not find testFile.ext")

  expect_error(getSamples(covFile,extFile=extFile),
               "n - number of samples - needs to be provided")

  set.seed("123")
  tmp <- getSamples(covFile,extFile=extFile,n=20)
  expect_is(tmp,"data.frame")
  expect_equal(nrow(tmp),21)

  ## Test that the output is the same as a saved reference
  expect_equal_to_reference(tmp,"test_output/getSamplesCovOutput")
})

test_that("getParamPositions works correctly", {

  bootFile  <- system.file("extdata","bs1_n100.dir/raw_results_run1.csv",package="PMXForest")
  rr_struct <- file.path(dirname(bootFile), "raw_results_structure")

  res <- getParamPositions(rr_struct)

  # ## Test that the output is the same as a saved reference
  expect_equal_to_reference(res,"test_output/getParamPositionBSOutput")
})


test_that("getSamples works correctly for bootstrap .csv input", {

  ## .csv file input
  bootFile <- system.file("extdata","bs1_n100.dir/raw_results_run1.csv",package="PMXForest")
  extFile  <- system.file("extdata","run1.ext",package="PMXForest")

  set.seed("123")
  tmp0 <- getSamples(bootFile,extFile=extFile)
  tmp1 <- getSamples(bootFile,extFile=extFile,n=200)

  ## Test that the output is the same as a saved reference
  expect_equal_to_reference(tmp0,"test_output/getSamplesBSOutput")
  expect_equal_to_reference(tmp1,"test_output/getSamplesBSN200Output")
})

test_that("getSamples works correctly for SIR .csv input", {

  ## .csv file input
  sirFile  <- system.file("extdata","sir_dir1/raw_results_run1.csv",package="PMXForest")
  extFile  <- system.file("extdata","run1.ext",package="PMXForest")

  set.seed("123")
  tmp0 <- getSamples(sirFile,extFile=extFile)

  ## Test that the output is the same as a saved reference
  expect_equal_to_reference(tmp0,"test_output/getSamplesSIROutput")
})


test_that("getSamples works correctly for data frame with bootstrap headings", {

  ## .csv file input
  bootFile <- system.file("extdata","bs1_n100.dir/raw_results_run1.csv",package="PMXForest")
  extFile  <- system.file("extdata","run1.ext",package="PMXForest")

  dft<-read.csv(bootFile)

  set.seed("123")
  tmp0 <- getSamples(subset(dft,ofv!=0))

  tmp1 <- getSamples(subset(dft,ofv!=0),indexvec = c(21:36,42,37:41))

  tmp2 <- getSamples(subset(dft,ofv!=0),indexvec = c(21:36,42,37:41),n=150)

  tmp3 <- getSamples(subset(dft,ofv!=0),indexvec = c(21:36,42,37:41),extFile = extFile)

  tmp4 <- getSamples(subset(dft,ofv!=0),indexvec = c(21:36,42,37:41),extFile = extFile,n=170)

  ## Test that the output is the same as a saved reference
  expect_equal_to_reference(tmp0,"test_output/dfSample0")
  expect_equal_to_reference(tmp1,"test_output/dfSample1")
  expect_equal_to_reference(tmp2,"test_output/dfSample2")
  expect_equal_to_reference(tmp3,"test_output/dfSample3")
  expect_equal_to_reference(tmp4,"test_output/dfSample4")
})


test_that("getSamples works correctly for a time-to-event model without SIGMA", {

  ## .csv file input
  runno<-"tte_weibull"
  bootFile <- system.file("extdata","tte","bootstrap_tte_weibull_n500",paste0("raw_results_",runno,".csv"),package="PMXForest")
  extFile <- system.file("extdata","tte",paste0(runno,".ext"),package="PMXForest")

  tmp <- getSamples(bootFile,extFile)

  ## Test that the output is the same as a saved reference
  expect_equal_to_reference(tmp,"test_output/dfSampleTTE")
})

test_that("getPlotVars returns correct values", {


  expect_equal(getPlotVars(plotRelative=FALSE,noVar=TRUE,reference="func"),c(point="POINT",q1="Q1",q2="Q2",ref="REFFUNC"))
  expect_equal(getPlotVars(plotRelative=FALSE,noVar=TRUE,reference="final"),c(point="POINT",q1="Q1",q2="Q2",ref="REFFINAL"))

  expect_error(getPlotVars(plotRelative=FALSE,noVar=FALSE,reference="func"))
  expect_error(getPlotVars(plotRelative=FALSE,noVar=FALSE,reference="final"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=TRUE,reference="func"),
               c(point="POINT_NOVAR_REL_REFFUNC",q1="Q1_NOVAR_REL_REFFUNC",q2="Q2_NOVAR_REL_REFFUNC",ref="REFFUNC"))

  expect_error(getPlotVars(plotRelative=TRUE,noVar=TRUE,reference="final"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=FALSE,reference="func"),
               c(point="POINT_REL_REFFUNC",q1="Q1_REL_REFFUNC",q2="Q2_REL_REFFUNC",ref="REFFUNC"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=FALSE,reference="final"),
               c(point="POINT_REL_REFFINAL",q1="Q1_REL_REFFINAL",q2="Q2_REL_REFFINAL",ref="REFFINAL"))

})
