library(testthat)

test_that("getPlotVars returns correct values", {

  expect_equal(getPlotVars(plotRelative=FALSE,noVar=TRUE,reference="func"),c(point="POINT",q1="Q1",q2="Q2",ref="REFFUNC"))
  expect_equal(getPlotVars(plotRelative=FALSE,noVar=TRUE,reference="final"),c(point="POINT",q1="Q1",q2="Q2",ref="REFTRUE"))

  expect_error(getPlotVars(plotRelative=FALSE,noVar=FALSE,reference="func"))
  expect_error(getPlotVars(plotRelative=FALSE,noVar=FALSE,reference="final"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=TRUE,reference="func"),
               c(point="POINT_NOVAR_REL_REFFUNC",q1="Q1_NOVAR_REL_REFFUNC",q2="Q2_NOVAR_REL_REFFUNC",ref="REFFUNC"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=TRUE,reference="final"),
               c(point="POINT_NOVAR_REL_REFTRUE",q1="Q1_NOVAR_REL_REFTRUE",q2="Q2_NOVAR_REL_REFTRUE",ref="REFTRUE"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=FALSE,reference="func"),
               c(point="POINT_REL_REFFUNC",q1="Q1_REL_REFFUNC",q2="Q2_REL_REFFUNC",ref="REFFUNC"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=FALSE,reference="final"),
               c(point="POINT_REL_REFTRUE",q1="Q1_REL_REFTRUE",q2="Q2_REL_REFTRUE",ref="REFTRUE"))

})

# test_that("setupForestPlotData works as expected", {
#
#   expect_is()
# })
