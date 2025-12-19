test_that("setupForestPlotData prepares data frame correctly", {
  # 1. Create mock dfres
  df_mock <- data.frame(
    PARAMETER = c("CL", "V"),
    GROUPNAME = c("Sex", "Weight"),
    COVNAME = c("Male", "70kg"),
    COVNUM = 1:2,
    COVEFF = c(TRUE, FALSE),
    REFROW = "NO",
    REFFUNC = 10,
    POINT_REL_REFFUNC = 1.234,
    Q1_REL_REFFUNC = 1.111,
    Q2_REL_REFFUNC = 1.357,
    stringsAsFactors = FALSE
  )

  # 2. Test Label Prefixing and Statistics (Lines 153-161, 203-215)
  res <- setupForestPlotData(
    df_mock,
    parameterLabelsPrefix = "Metric: ",
    sigdigits = 2
  )

  expect_s3_class(res, "data.frame")
  # Check Prefixing
  expect_equal(as.character(res$PARAMETERLABEL[1]), "Metric: CL")
  # Check Statistic String Padding
  # 1.234 to 2 sig digits is 1.2
  expect_match(res$STATISTIC[1], "1.2 \\[1.1-1.4\\]")

  # 3. Test Significance Filtering (Lines 181-183)
  res_sig <- setupForestPlotData(df_mock, onlySignificant = TRUE)
  # Should only contain 'Sex' because 'Weight' COVEFF is FALSE
  expect_equal(nrow(res_sig), 1)
  expect_equal(res_sig$GROUPNAME, "Sex")

  # 4. Test Group Name Custom Labels (Lines 166-173)
  res_labels <- setupForestPlotData(df_mock, groupNameLabels = c("Gender", "Body Mass"))
  expect_equal(as.character(res_labels$GROUPNAMELABEL[1]), "Gender")
  expect_equal(as.character(res_labels$GROUPNAMELABEL[2]), "Body Mass")
})

test_that("setupForestPlotData validation errors", {
  df_mock <- data.frame(PARAMETER = "CL", GROUPNAME = "Sex", stringsAsFactors = FALSE)

  # Trigger Parameter Label Error (Line 131)
  expect_error(setupForestPlotData(df_mock, parameters = "CL", parameterLabels = c("L1", "L2")),
               "number of parameter labels must either be the same")

  # Trigger Group Name Label Error (Line 138)
  expect_error(setupForestPlotData(df_mock, groupNameLabels = c("G1", "G2")),
               "number of group name labels must either be the same")
})

test_that("setupForestPlotData hits remaining logical branches", {
  # Setup deterministic mock data
  df_mock <- data.frame(
    PARAMETER = c("CL", "V"),
    GROUPNAME = c("Sex", "Weight"),
    COVNAME   = c("Male", "70kg"),
    COVNUM    = 1:2,
    COVEFF    = c(TRUE, FALSE),
    REFROW    = "NO",
    REFFUNC   = 10,
    POINT_REL_REFFUNC = 1.2,
    Q1_REL_REFFUNC    = 1.1,
    Q2_REL_REFFUNC    = 1.3,
    stringsAsFactors  = FALSE
  )

  # 1. Trigger Custom Parameter and Statistics Labels (Lines 77-78, 87)
  res_custom_labels <- setupForestPlotData(
    df_mock,
    parameterLabels  = c("Clearance", "Volume"),
    statisticsLabels = "Stats: "
  )
  expect_equal(as.character(res_custom_labels$PARAMETERLABEL[1]), "Clearance")
  expect_equal(as.character(res_custom_labels$STATISTICSLABEL[1]), "Stats: Clearance")

  # 2. Trigger Manual Group Name Labels (Line 109)
  # Providing a vector of labels that matches the row count instead of unique group count
  res_group_labels <- setupForestPlotData(
    df_mock,
    groupNameLabels = c("LabelRow1", "LabelRow2")
  )
  expect_equal(as.character(res_group_labels$GROUPNAMELABEL[1]), "LabelRow1")

  # 3. Trigger Significance Overrides (Line 125)
  # Using setSignEff to manually flip COVEFF for the Weight group
  # Note: Requires setCOVEFF helper function to be available
  lsSignEff <- list(c("V", "Weight"))
  res_override <- setupForestPlotData(df_mock, setSignEff = lsSignEff)

  # Check if the override flipped COVEFF to TRUE for the second row
  expect_true(res_override$COVEFF[2])
})
