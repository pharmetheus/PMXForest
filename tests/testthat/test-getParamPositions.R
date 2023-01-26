test_that("getParamPositions works correctly", {

  bootFile  <- system.file("extdata","SimVal/bs7.dir/raw_results_run7bs.csv",package="PMXForest")
  rr_struct <- file.path(dirname(bootFile), "raw_results_structure")

  res <- getParamPositions(rr_struct)

  # ## Test that the output is the same as a saved reference
  expect_snapshot(res)
})
