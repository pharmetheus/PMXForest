context("Sampling for uncertainty estimation")
library(testthat)
library(PMXForest)


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

  ## Test iSampleIndex
  dfParameters <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = 10))
  expect_equal(nrow(dfParameters),10)

  dfParameters <- mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = 1)
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
  expect_equal(nrow(tmp),20)

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

