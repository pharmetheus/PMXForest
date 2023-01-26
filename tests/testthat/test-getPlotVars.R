test_that("getPlotVars returns correct values", {

  expect_error(getPlotVars(reference="tmp"))

  expect_equal(getPlotVars(plotRelative=FALSE,noVar=TRUE,reference="func"),c(REF="REFFUNC",point="POINT",q1="Q1",q2="Q2"))
  expect_equal(getPlotVars(plotRelative=FALSE,noVar=TRUE,reference="final"),c(REF="REFFINAL",point="POINT",q1="Q1",q2="Q2"))

  expect_error(getPlotVars(plotRelative=FALSE,noVar=FALSE,reference="func"))
  expect_error(getPlotVars(plotRelative=FALSE,noVar=FALSE,reference="final"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=TRUE,reference="func"),
               c(REF="REFFUNC",point="POINT_NOVAR_REL_REFFUNC",q1="Q1_NOVAR_REL_REFFUNC",q2="Q2_NOVAR_REL_REFFUNC"))

  expect_error(getPlotVars(plotRelative=TRUE,noVar=TRUE,reference="final"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=FALSE,reference="func"),
               c(REF="REFFUNC",point="POINT_REL_REFFUNC",q1="Q1_REL_REFFUNC",q2="Q2_REL_REFFUNC"))

  expect_equal(getPlotVars(plotRelative=TRUE,noVar=FALSE,reference="final"),
               c(REF="REFFINAL",point="POINT_REL_REFFINAL",q1="Q1_REL_REFFINAL",q2="Q2_REL_REFFINAL"))

})




