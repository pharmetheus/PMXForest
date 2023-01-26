test_that("Ext files are read properly", {
  extData <- getExt(system.file("extdata","SimVal/run7.ext",package="PMXForest"))
  expect_type(extData,"list")
})
