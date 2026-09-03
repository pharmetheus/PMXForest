modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")

## Write a NONMEM-format table: a "TABLE NO." banner, a header, then the rows.
writeNmTable <- function(df, path) {
  writeLines(c(
    "TABLE NO.  1",
    paste0(" ", paste(formatC(names(df), width = 11, flag = "-"), collapse = " ")),
    apply(df, 1, function(r)
      paste0(" ", paste(formatC(as.numeric(r), format = "E", digits = 4, width = 11),
                        collapse = " ")))
  ), path)
  path
}

## A table built from the generated function itself, so a correct
## implementation must reproduce it exactly. Verification is being tested here,
## not the translation - that is covered in test-createParamFunction.R.
makeTable <- function(out, thetas, n = 20, withEta = TRUE, withCovs = TRUE) {
  set.seed(42)
  covs <- data.frame(
    ID    = seq_len(n),
    WT    = seq(55, 120, length.out = n),
    SEX   = rep(c(1, 2), length.out = n),
    FOOD  = rep(c(0, 1), length.out = n),
    FORM  = rep(c(0, 1), length.out = n),
    GENO1 = rep(c(0, 1, 0, 0), length.out = n),
    GENO3 = rep(c(0, 0, 1, 0), length.out = n),
    GENO4 = rep(c(0, 0, 0, 1), length.out = n)
  )
  fun <- eval(parse(text = out$code))
  # rbind rather than t(vapply()), which would transpose for a single parameter
  tv <- do.call(rbind, lapply(seq_len(n), function(i)
    unlist(fun(thetas, covs[i, , drop = FALSE]))))

  tab <- covs
  etas <- stats::rnorm(n, sd = 0.2)
  for (p in out$functionListName) {
    idx <- out$etaMap[[p]]
    if (withEta && !is.null(idx)) {
      # NONMEM tables individual values: TV * exp(eta)
      tab[[p]] <- tv[, p] * exp(etas)
      tab[[paste0("ETA", idx)]] <- etas
    } else {
      tab[[p]] <- tv[, p]
    }
  }
  if (!withCovs) tab$WT <- NULL
  tab
}

run7Thetas <- function() {
  ext <- getExt(system.file("extdata", "SimVal/run7.ext", package = "PMXForest"))
  as.numeric(ext[ext$ITERATION == -1000000000, 2:15])
}

test_that("a faithful function passes against the table it produced", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE)
  )
  thetas <- run7Thetas()
  f   <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(makeTable(out, thetas), f)

  v <- verifyParamFunction(out, f, thetas, quiet = TRUE)

  # scalar TRUE/FALSE, usable directly in an if
  expect_length(as.logical(v), 1L)
  expect_true(as.logical(v))
  expect_true(if (v) TRUE else FALSE)

  d <- attr(v, "checks")
  expect_s3_class(d, "data.frame")
  expect_named(d, c("PARAMETER", "TABLECOLUMN", "N", "MAXABSDIFF",
                    "MAXRELDIFF", "PASS"))
  expect_equal(d$PARAMETER, c("CL", "V"))
  expect_true(all(d$PASS))
  # the individual values were reconciled by dividing out the exponential IIV
  expect_match(d$TABLECOLUMN[1], "CL / exp\\(ETA3\\)")
  expect_true(all(d$MAXRELDIFF < 1e-4))
})

test_that("a wrong function fails", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE)
  )
  thetas <- run7Thetas()
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(makeTable(out, thetas), f)

  # A deliberately mistranslated function: the WT exponent uses the wrong theta.
  broken <- eval(parse(text = sub("thetas\\[2\\]", "thetas[3]",
                                  paste(out$code, collapse = "\n"))))
  v <- verifyParamFunction(out, f, thetas, fun = broken, quiet = TRUE)
  expect_false(as.logical(v))
  d <- attr(v, "checks")
  expect_false(d$PASS[d$PARAMETER == "CL"])
  expect_true(d$MAXRELDIFF[d$PARAMETER == "CL"] > 1e-4)
})

test_that("a TV-prefixed column is used directly, without needing an ETA", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  thetas <- run7Thetas()
  tab <- makeTable(out, thetas, withEta = FALSE)
  names(tab)[names(tab) == "CL"] <- "TVCL"

  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(tab, f)
  v <- verifyParamFunction(out, f, thetas, quiet = TRUE)
  d <- attr(v, "checks")
  expect_equal(d$TABLECOLUMN, "TVCL")
  expect_true(as.logical(v))
  expect_true(d$PASS)
})

test_that("a missing covariate warns loudly and voids the result", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  thetas <- run7Thetas()
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(makeTable(out, thetas, withCovs = FALSE), f)

  expect_warning(v <- verifyParamFunction(out, f, thetas, quiet = TRUE),
                 "not columns of")
  # PASS is NA rather than TRUE/FALSE: the comparison is not a valid check,
  # and an unverifiable parameter makes the scalar FALSE.
  expect_false(as.logical(v))
  expect_true(is.na(attr(v, "checks")$PASS))
})

test_that("a parameter absent from the table warns and is reported unverified", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("CL", "MAT"), quiet = TRUE)
  )
  thetas <- run7Thetas()
  tab <- makeTable(out, thetas)
  tab$MAT <- NULL
  tab$ETA5 <- NULL

  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(tab, f)
  expect_warning(v <- verifyParamFunction(out, f, thetas, quiet = TRUE),
                 "not a column of")
  d <- attr(v, "checks")
  expect_true(d$PASS[d$PARAMETER == "CL"])
  expect_true(is.na(d$PASS[d$PARAMETER == "MAT"]))
  expect_equal(d$N[d$PARAMETER == "MAT"], 0L)
  # one parameter could not be verified -> overall FALSE
  expect_false(as.logical(v))
})

test_that("an individual column with no ETA column warns", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  thetas <- run7Thetas()
  tab <- makeTable(out, thetas)
  tab$ETA3 <- NULL          # CL is tabled but cannot be reduced to typical

  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(tab, f)
  expect_warning(v <- verifyParamFunction(out, f, thetas, quiet = TRUE),
                 "Cannot recover typical values")
  expect_true(is.na(attr(v, "checks")$PASS))
  expect_false(as.logical(v))
})

test_that("rows are deduplicated on the covariate combination", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  thetas <- run7Thetas()
  tab <- makeTable(out, thetas, n = 5)
  tab <- tab[rep(seq_len(nrow(tab)), each = 4), ]   # 4 records per subject

  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(tab, f)
  v <- verifyParamFunction(out, f, thetas, quiet = TRUE)
  expect_equal(attr(v, "checks")$N, 5L)
})

test_that("the detail attribute carries the per-row comparison", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  thetas <- run7Thetas()
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(makeTable(out, thetas, n = 8), f)

  v      <- verifyParamFunction(out, f, thetas, quiet = TRUE)
  detail <- attr(attr(v, "checks"), "detail")
  expect_named(detail, "CL")
  expect_named(detail$CL, c("generated", "table", "absdiff", "reldiff"))
  expect_equal(nrow(detail$CL), 8)
})

test_that("printing shows the per-parameter table", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  thetas <- run7Thetas()
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(makeTable(out, thetas), f)

  v <- verifyParamFunction(out, f, thetas, quiet = TRUE)
  expect_output(print(v), "PASS - verifyParamFunction: 1/1")
  expect_output(print(v), "CL")
})

test_that("quiet = FALSE reports each parameter", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  thetas <- run7Thetas()
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(makeTable(out, thetas), f)
  expect_message(verifyParamFunction(out, f, thetas), "CL: pass")
})

test_that("bad input is rejected", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  expect_error(verifyParamFunction(list(), "nowhere.tab", 1),
               "returned by createParamFunction")
  expect_error(verifyParamFunction(out, "does-not-exist.tab", 1),
               "Table file not found")
})
