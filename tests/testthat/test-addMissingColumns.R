test_that("addMissingColumns transforms parameters correctly", {
  # 1. Setup mock .ext template
  # Physical Positions:
  # 1: ITERATION
  # 2: THETA1
  # 3: OMEGA.3.1. (Value 0, target for detection)
  # 4: OMEGA.3.2. (Value 0, target for detection)
  # 5: OMEGA.5.1. (Non-zero)
  # 6: OBJ
  df_ext <- data.frame(
    ITERATION  = -1e9,
    THETA1     = 1.0,
    OMEGA.3.1. = 0.0,
    OMEGA.3.2. = 0.0,
    OMEGA.5.1. = 0.5,
    OBJ        = -5000,
    stringsAsFactors = FALSE
  )

  # 2. Setup df_params
  # Template has 4 parameter columns (THETA1, O31, O32, O51).
  # Two are 'missing' (O31, O32) because they are zero in the template.
  # Therefore, df_params must have exactly 2 columns.
  df_params <- data.frame(
    V1 = 1.1,
    V2 = 0.6,
    stringsAsFactors = FALSE
  )

  # Case A: Automatic zerosindex detection (Lines 23-32)
  # This triggers the nested loops to find indices 3 and 4.
  res_auto <- addMissingColumns(df_params, df_ext, zerosindex = NULL)

  expect_s3_class(res_auto, "data.frame")
  expect_equal(ncol(res_auto), 4) # Recovers the two missing columns
  expect_equal(res_auto$OMEGA.3.1.[1], 0.0)
  expect_equal(res_auto$THETA1[1], 1.1)

  # Case B: Manual zerosindex (Line 36)
  # Directly provide the indices for the zero columns (3 and 4)
  res_manual <- addMissingColumns(df_params, df_ext, zerosindex = c(3, 4))
  expect_equal(res_manual$OMEGA.3.2.[1], 0.0)
})

test_that("addMissingColumns handles no fixed parameters", {
  # Trigger the 'else' branch (Line 46)
  df_ext <- data.frame(
    ITERATION = -1e9,
    THETA1 = 1.0,
    OBJ = -5000,
    stringsAsFactors = FALSE
  )
  df_params <- data.frame(V1 = 1.1, stringsAsFactors = FALSE)

  res_none <- addMissingColumns(df_params, df_ext, zerosindex = NULL)
  expect_named(res_none, "THETA1")
})
