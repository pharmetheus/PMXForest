test_that("Ext files are read properly", {
  extData <- getExt(system.file("extdata","SimVal/run7.ext",package="PMXForest"))
  expect_s3_class(extData, "data.frame")
})

test_that("getExt handles all logic branches", {
  # Use a temporary directory that is automatically cleaned up
  temp_dir <- withr::local_tempdir()

  create_file <- function(path, content) {
    full_path <- file.path(temp_dir, path)
    writeLines(text = content, con = full_path)
    return(full_path)
  }

  # --- Test Fixtures ---
  file_1_table <- create_file("one_table.ext", c("TABLE NO. 1", "ITERATION TVAL", "-1000000000 1.1"))
  file_2_tables <- create_file("two_tables.ext", c("TABLE NO. 1", "ITERATION TVAL", "-1000000000 1.1", "TABLE NO. 2", "ITERATION TVAL", "-1000000000 2.2"))
  file_3_tables <- create_file("three_tables.ext", c("TABLE NO. 1", "ITERATION TVAL", "-1000000000 1.1", "TABLE NO. 2", "ITERATION TVAL", "-1000000000 2.2", "TABLE NO. 3", "ITERATION TVAL", "-1000000000 3.3"))
  file_4_tables <- create_file("four_tables.ext", c("TABLE NO. 1", "A B", "1 1", "TABLE NO. 2", "C D", "2 2", "TABLE NO. 3", "E F", "3 3", "TABLE NO. 4", "ITERATION TVAL", "-1000000000 4.4"))
  file_no_table <- create_file("no_table.ext", c("HEADER LINE", "ITERATION TVAL", "1 1.1"))

  # --- Test Cases ---
  # Test the default behavior (set=NULL, should get the last table)
  expect_equal(getExt(file_1_table)$TVAL, 1.1)
  expect_equal(getExt(file_4_tables)$TVAL, 4.4)

  # Test each specific logic branch from the function
  expect_equal(getExt(file_2_tables, set = 1)$TVAL, 1.1) # set=1, tables > 1
  expect_equal(getExt(file_2_tables, set = 2)$TVAL, 2.2) # set=2, tables = 2
  expect_equal(getExt(file_3_tables, set = 2)$TVAL, 2.2) # set=2, tables = 3
  expect_equal(getExt(file_3_tables, set = 3)$TVAL, 3.3) # set=3, tables = 3
  expect_equal(getExt(file_4_tables, set = 4)$TVAL, 4.4) # set=4, tables = 4

  # Test error conditions
  expect_error(getExt(file_1_table, set = 2)) # Requesting a table that doesn't exist
  expect_error(getExt(file_no_table))         # File contains no "TABLE" headers
})
