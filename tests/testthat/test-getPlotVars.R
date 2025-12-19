test_that("getPlotVars returns correct values and handles errors", {

  # --- Error Handling ---
  # Test invalid reference
  expect_error(getPlotVars(reference = "tmp"),
               "reference needs to be either func or final.")

  # Test invalid combinations
  expect_error(getPlotVars(plotRelative = FALSE, noVar = FALSE, reference = "func"),
               "combination of plotRelative=FALSE and noVar=FALSE is not possible")

  expect_error(getPlotVars(plotRelative = TRUE, noVar = TRUE, reference = "final"),
               "combination of plotRelative=TRUE, noVar=TRUE and reference=final is not possible")

  # --- Valid Combinations (R-version independent) ---
  # Non-relative with func
  expect_equal(getPlotVars(plotRelative = FALSE, noVar = TRUE, reference = "func"),
               c(REF = "REFFUNC", point = "POINT", q1 = "Q1", q2 = "Q2"))

  # Non-relative with final
  expect_equal(getPlotVars(plotRelative = FALSE, noVar = TRUE, reference = "final"),
               c(REF = "REFFINAL", point = "POINT", q1 = "Q1", q2 = "Q2"))

  # Relative with noVar and func
  expect_equal(getPlotVars(plotRelative = TRUE, noVar = TRUE, reference = "func"),
               c(REF = "REFFUNC", point = "POINT_NOVAR_REL_REFFUNC",
                 q1 = "Q1_NOVAR_REL_REFFUNC", q2 = "Q2_NOVAR_REL_REFFUNC"))

  # Relative with uncertainty (func)
  expect_equal(getPlotVars(plotRelative = TRUE, noVar = FALSE, reference = "func"),
               c(REF = "REFFUNC", point = "POINT_REL_REFFUNC",
                 q1 = "Q1_REL_REFFUNC", q2 = "Q2_REL_REFFUNC"))

  # Relative with uncertainty (final)
  expect_equal(getPlotVars(plotRelative = TRUE, noVar = FALSE, reference = "final"),
               c(REF = "REFFINAL", point = "POINT_REL_REFFINAL",
                 q1 = "Q1_REL_REFFINAL", q2 = "Q2_REL_REFFINAL"))
})



