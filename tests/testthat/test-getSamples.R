test_that("getSamples works correctly for .cov input", {

  ## .cov file input
  covFile <- system.file("extdata","SimVal/run7.cov",package="PMXForest")
  extFile <- system.file("extdata","SimVal/run7.ext",package="PMXForest")

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
  expect_type(tmp,"list")
  expect_equal(nrow(tmp),21)

  ## Test that the output is the same as a saved reference
  expect_snapshot(tmp)
})

test_that("getSamples works correctly for bootstrap .csv input", {

  ## .csv file input
  bootFile <- system.file("extdata","SimVal/bs7.dir/raw_results_run7bs.csv",package="PMXForest")
  extFile  <- system.file("extdata","SimVal/run7.ext",package="PMXForest")

  set.seed("123")
  expect_error(getSamples(bootFile))
  tmp0 <- getSamples(bootFile,extFile=extFile)
  tmp1 <- getSamples(bootFile,extFile=extFile,n=200)

  tmp2<- getSamples(bootFile,extFile=extFile,indexvec = c(21:34,40,35:39))
  expect_snapshot(tmp2)

  ## Test that the output is the same as a saved reference
  expect_snapshot(tmp0)
  expect_snapshot(tmp1)
})

test_that("getSamples works correctly for SIR .csv input", {

  ## .csv file input
  sirFile  <- system.file("extdata","SimVal/sir7.dir/raw_results_run7.csv",package="PMXForest")
  extFile  <- system.file("extdata","SimVal/run7.ext",package="PMXForest")

  set.seed("123")
  tmp0 <- getSamples(sirFile,extFile=extFile)

  ## Test that the output is the same as a saved reference
  expect_snapshot(tmp0)
})


test_that("getSamples works correctly for data frame with bootstrap headings", {

  ## .csv file input
  bootFile <- system.file("extdata","SimVal/bs7.dir/raw_results_run7bs.csv",package="PMXForest")
  extFile  <- system.file("extdata","SimVal/run7.ext",package="PMXForest")

  dft <- read.csv(bootFile,stringsAsFactors = FALSE)

  set.seed("123")
  tmp0 <- getSamples(subset(dft,ofv!=0))

  ## Check that ext-file is provided
  tmp1 <- getSamples(subset(dft,ofv!=0),indexvec = c(21:34,40,35:39))

  tmp2 <- getSamples(subset(dft,ofv!=0),indexvec = c(21:34,40,35:39),n=150)

  tmp3 <- getSamples(subset(dft,ofv!=0),indexvec = c(21:34,40,35:39),extFile = extFile)

  tmp4 <- getSamples(subset(dft,ofv!=0),indexvec = c(21:34,40,35:39),extFile = extFile,n=170)

  ## Test that the output is the same as a saved reference
  expect_snapshot(tmp0)
  expect_snapshot(tmp1)
  expect_snapshot(tmp2)
  expect_snapshot(tmp3)
  expect_snapshot(tmp4)
})


test_that("getSamples works correctly for a time-to-event model without SIGMA", {

  ## .csv file input
  runno<-"tte_weibull"
  bootFile <- system.file("extdata","tte","bootstrap_tte_weibull_n500",paste0("raw_results_",runno,".csv"),package="PMXForest")
  extFile <- system.file("extdata","tte",paste0(runno,".ext"),package="PMXForest")

  tmp <- getSamples(bootFile,extFile)

  ## Test that the output is the same as a saved reference
  expect_snapshot(tmp)
})

