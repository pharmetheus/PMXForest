test_that("addStamp adds a caption and returns a ggplot", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point()
  s <- addStamp(p)

  expect_s3_class(s, "ggplot")
  expect_match(s$labels$caption, "^.+\nCreated:")
  # the original plot is not modified
  expect_null(p$labels$caption)
})

test_that("addStamp records the working directory and the time", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point()
  cap <- addStamp(p)$labels$caption

  expect_match(cap, basename(getwd()), fixed = TRUE)
  # a parseable timestamp on the second line
  stamped <- sub("^.*\nCreated:", "", cap)
  expect_false(is.na(as.POSIXct(stamped)))
})

test_that("addStamp honours an explicit source and size", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point()
  cap <- addStamp(p, source = "somewhere/else.Rmd/fig-1")$labels$caption

  expect_match(cap, "^somewhere/else\\.Rmd/fig-1\nCreated:")
  expect_no_match(cap, basename(getwd()), fixed = TRUE)

  s <- addStamp(p, size = 11)
  expect_equal(s$theme$plot.caption$size, 11)
})

test_that("addStamp rejects a non-ggplot", {
  expect_error(addStamp(mtcars), "must be a ggplot object")
  expect_error(addStamp("not a plot"), "must be a ggplot object")
})

test_that("the caption omits unknown segments rather than leaving empty ones", {
  # outside a knit session there is no input file and no chunk label, so the
  # source is the working directory alone - not "dir//" with empty segments
  txt <- PMXForest:::stampText()
  expect_equal(sub("\nCreated:.*$", "", txt), basename(getwd()))
  expect_no_match(txt, "//", fixed = TRUE)
})
