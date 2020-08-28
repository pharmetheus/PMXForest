context("Utility functions")
library(testthat)
library(PMXForest)


test_that("Ext files are read properly", {
  extData <- getExt(system.file("extdata","run1.ext",package="PMXForest"))
  expect_is(extData,"data.frame")
})
