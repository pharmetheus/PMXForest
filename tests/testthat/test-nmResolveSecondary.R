test_that("NULL / empty give an empty list", {
  expect_identical(nmResolveSecondary(NULL), list())
  expect_identical(nmResolveSecondary(list()), list())
})

test_that("a snippet is classified and split into lines", {
  r <- nmResolveSecondary(list(AUC = "df$DOSE / CL"), quiet = TRUE)
  expect_named(r, "AUC")
  expect_equal(r$AUC$name, "AUC")
  expect_true(is.na(r$AUC$src))
  expect_equal(r$AUC$lines, "df$DOSE / CL")

  r2 <- nmResolveSecondary(list(X = "a <- CL / V\na * 2"), quiet = TRUE)
  expect_equal(r2$X$lines, c("a <- CL / V", "a * 2"))
})

test_that("an existing file is read; a missing .R path is an error", {
  f <- withr::local_tempfile(fileext = ".R")
  writeLines(c("k <- CL / V", "k"), f)
  r <- nmResolveSecondary(list(K = f), quiet = TRUE)
  expect_equal(r$K$src, f)
  expect_equal(r$K$lines, c("k <- CL / V", "k"))

  expect_error(nmResolveSecondary(list(K = "no-such-file.R")),
               "looks like a file path but does not exist")
})

test_that("a snippet that happens to have no .R extension is not mistaken for a path", {
  # 'CL' is not a file and has no .R extension -> snippet, no error
  r <- nmResolveSecondary(list(P = "CL / V"), quiet = TRUE)
  expect_true(is.na(r$P$src))
})

test_that("names must be present, unique and valid", {
  expect_error(nmResolveSecondary(list("CL / V")), "named list")
  expect_error(nmResolveSecondary(list(A = "1", A = "2")), "duplicate")
  expect_error(nmResolveSecondary(list(`a b` = "1")), "not valid R names")
})

test_that("a non-string entry is rejected", {
  expect_error(nmResolveSecondary(list(A = 1)), "single string")
  expect_error(nmResolveSecondary(list(A = c("a", "b"))), "single string")
})

test_that("a syntax error names the offending entry", {
  expect_error(nmResolveSecondary(list(BAD = "CL / ")), "secondary 'BAD'")
  f <- withr::local_tempfile(fileext = ".R")
  writeLines("function(", f)
  expect_error(nmResolveSecondary(list(BAD = f)), "secondary 'BAD'")
})

test_that("quiet = FALSE reports how each entry was read", {
  f <- withr::local_tempfile(fileext = ".R")
  writeLines("CL", f)
  expect_message(nmResolveSecondary(list(A = "CL / V")), "secondary A: inline snippet")
  expect_message(nmResolveSecondary(list(B = f)), "inlined from ")
})
