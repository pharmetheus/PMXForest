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
  expect_error(
    setupForestPlotData(df_mock, parameters = "CL", parameterLabels = c("L1", "L2")),
    "number of parameter labels must either be the same"
  )

  # Trigger Group Name Label Error (Line 138)
  expect_error(
    setupForestPlotData(df_mock, groupNameLabels = c("G1", "G2")),
    "number of group name labels must either be the same"
  )
})

test_that("setupForestPlotData hits remaining logical branches", {
  # Setup deterministic mock data
  df_mock <- data.frame(
    PARAMETER = c("CL", "V"),
    GROUPNAME = c("Sex", "Weight"),
    COVNAME = c("Male", "70kg"),
    COVNUM = 1:2,
    COVEFF = c(TRUE, FALSE),
    REFROW = "NO",
    REFFUNC = 10,
    POINT_REL_REFFUNC = 1.2,
    Q1_REL_REFFUNC = 1.1,
    Q2_REL_REFFUNC = 1.3,
    stringsAsFactors = FALSE
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

test_that("forestPlot's default statisticsLabel separates itself from the parameter", {
  # The label is prepended verbatim (see the "Stats: " case above), so the
  # separating space has to be part of the default. Without it the facet strip
  # read "Statistics:CL (L/h)" - which is what the README figure showed for as
  # long as that default lacked the space.
  expect_equal(eval(formals(forestPlot)$statisticsLabel), "Statistics: ")

  df_mock <- data.frame(
    PARAMETER = c("CL", "V"), GROUPNAME = c("Sex", "Weight"),
    COVNAME = c("Male", "70kg"), COVNUM = 1:2, COVEFF = c(TRUE, FALSE),
    REFROW = "NO", REFFUNC = 10, POINT_REL_REFFUNC = 1.234,
    Q1_REL_REFFUNC = 1.111, Q2_REL_REFFUNC = 1.357, stringsAsFactors = FALSE
  )
  res <- setupForestPlotData(
    df_mock,
    statisticsLabels = eval(formals(forestPlot)$statisticsLabel)
  )
  expect_equal(as.character(res$STATISTICSLABEL[1]), "Statistics: CL")
  expect_false(grepl("Statistics:[^ ]", as.character(res$STATISTICSLABEL[1])))
})

test_that("setupForestPlotData accepts groupNameLabels as a per-row vector", {
  df_mock <- data.frame(
    PARAMETER = c("CL", "CL", "CL"),
    GROUPNAME = c("Sex", "Sex", "Weight"), # 2 unique groups, 3 rows
    COVNAME = c("Male", "Female", "70kg"),
    COVNUM = 1:3,
    COVEFF = TRUE,
    REFROW = "NO",
    REFFUNC = 10,
    POINT_REL_REFFUNC = 1.2, Q1_REL_REFFUNC = 1.1, Q2_REL_REFFUNC = 1.3,
    stringsAsFactors = FALSE
  )

  res <- setupForestPlotData(df_mock, groupNameLabels = c("A", "B", "C"))
  expect_equal(as.character(res$GROUPNAMELABEL), c("A", "B", "C"))
})

# --- statistics-table number format: sigdigits vs decimals ---

test_that("setupForestPlotData chooses decimals on the relative scale, sigdigits on the absolute", {
  df_mock <- data.frame(
    PARAMETER = "CL", GROUPNAME = "Weight", COVNAME = "115 kg",
    COVNUM = 1, COVEFF = TRUE, REFROW = "NO", REFFUNC = 10, REFFINAL = 10,
    POINT = 16.401, Q1 = 15.62, Q2 = 17.284,
    POINT_NOVAR_REL_REFFUNC = 1.234, Q1_NOVAR_REL_REFFUNC = 1.151,
    Q2_NOVAR_REL_REFFUNC = 1.357,
    stringsAsFactors = FALSE
  )

  # relative, default -> 2 decimals (1.234 keeps the second decimal, not "1.2")
  rel <- setupForestPlotData(df_mock, plotRelative = TRUE, noVar = TRUE)
  expect_equal(trimws(rel$STATISTIC[1]), "1.23 [1.15-1.36]")

  # absolute, default -> 2 significant digits (unchanged behaviour)
  abs <- setupForestPlotData(df_mock,
    plotRelative = FALSE, noVar = TRUE,
    reference = "func"
  )
  expect_equal(trimws(abs$STATISTIC[1]), "16 [16-17]")
})

test_that("setupForestPlotData honours explicit sigdigits / decimals on either scale", {
  df_mock <- data.frame(
    PARAMETER = "CL", GROUPNAME = "Weight", COVNAME = "115 kg",
    COVNUM = 1, COVEFF = TRUE, REFROW = "NO", REFFUNC = 10, REFFINAL = 10,
    POINT = 16.401, Q1 = 15.62, Q2 = 17.284,
    POINT_NOVAR_REL_REFFUNC = 1.234, Q1_NOVAR_REL_REFFUNC = 1.151,
    Q2_NOVAR_REL_REFFUNC = 1.357,
    stringsAsFactors = FALSE
  )

  # explicit sigdigits on the relative scale
  expect_equal(
    trimws(setupForestPlotData(df_mock,
      plotRelative = TRUE, noVar = TRUE,
      sigdigits = 3
    )$STATISTIC[1]),
    "1.23 [1.15-1.36]"
  )
  # explicit decimals on the absolute scale
  expect_equal(
    trimws(setupForestPlotData(df_mock,
      plotRelative = FALSE, noVar = TRUE,
      reference = "func", decimals = 1
    )$STATISTIC[1]),
    "16.4 [15.6-17.3]"
  )
  # both is an error
  expect_error(
    setupForestPlotData(df_mock, sigdigits = 2, decimals = 2),
    "either .sigdigits. or .decimals."
  )
})

test_that("signifPad rounds to significant digits half-up and pads trailing zeros", {
  sp <- PMXForest:::signifPad

  expect_equal(
    sp(c(0.976, 1.234, 2.244, 12.3, 0.08), digits = 3),
    c("0.976", "1.23", "2.24", "12.3", "0.0800")
  )
  expect_equal(sp(c(1.2, 1.234, 1.15), digits = 2), c("1.2", "1.2", "1.2"))
  expect_equal(sp(c(16.4, 15.62, 17.284), digits = 2), c("16", "16", "17"))
  expect_equal(sp(-0.5, digits = 2), "-0.50")
  expect_equal(sp(100, digits = 2), "100") # no bare trailing "."
  expect_true(is.na(sp(NA_real_, digits = 2)))
})
