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

  expect_error(
    nmResolveSecondary(list(K = "no-such-file.R")),
    "looks like a file path but does not exist"
  )
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

test_that("a non-string, non-list entry is rejected", {
  expect_error(nmResolveSecondary(list(A = 1)), "single string")
  expect_error(nmResolveSecondary(list(A = c("a", "b"))), "single string")
})

test_that("a list entry: `source` plus constants become `name <- value` lines", {
  r <- nmResolveSecondary(
    list(CMAX = list(
      source = "cmax(dose, tau)", dose = 100, tau = 12,
      label = "qd", flag = TRUE, times = c(0, 24, 48)
    )),
    quiet = TRUE
  )
  expect_equal(r$CMAX$name, "CMAX")
  expect_true(is.na(r$CMAX$src)) # source is a snippet here
  expect_equal(r$CMAX$lines, "cmax(dose, tau)")
  expect_equal(
    r$CMAX$consts,
    c(
      "dose <- 100", "tau <- 12", "label <- \"qd\"",
      "flag <- TRUE", "times <- c(0, 24, 48)"
    )
  )
})

test_that("a list entry with a file `source` is read, constants kept", {
  f <- withr::local_tempfile(fileext = ".R")
  writeLines(c("dose / CL"), f)
  r <- nmResolveSecondary(list(AUC = list(source = f, dose = 50)), quiet = TRUE)
  expect_equal(r$AUC$src, f)
  expect_equal(r$AUC$lines, "dose / CL")
  expect_equal(r$AUC$consts, "dose <- 50")
})

test_that("a bare string entry has empty consts", {
  r <- nmResolveSecondary(list(AUC = "df$DOSE / CL"), quiet = TRUE)
  expect_identical(r$AUC$consts, character(0))
})

test_that("a list entry is validated", {
  expect_error(
    nmResolveSecondary(list(A = list(dose = 100))),
    "needs a single-string `source`"
  )
  expect_error(
    nmResolveSecondary(list(A = list(source = "x", 100))),
    "must be .*named"
  )
  expect_error(
    nmResolveSecondary(list(A = list(source = "x", `a b` = 1))),
    "not valid R names"
  )
  expect_error(
    nmResolveSecondary(list(A = list(source = "x", d = list(1)))),
    "must be an atomic value"
  )
})

test_that("quiet = FALSE notes how many constants an entry carries", {
  expect_message(
    nmResolveSecondary(list(A = list(source = "dose / CL", dose = 100))),
    "\\(\\+1 constant"
  )
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
