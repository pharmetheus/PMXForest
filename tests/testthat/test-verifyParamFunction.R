modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")

## Write a NONMEM-format table: a "TABLE NO." banner, a header, then the rows.
writeNmTable <- function(df, path) {
  writeLines(c(
    "TABLE NO.  1",
    paste0(" ", paste(formatC(names(df), width = 11, flag = "-"), collapse = " ")),
    apply(df, 1, function(r) {
      paste0(" ", paste(formatC(as.numeric(r), format = "E", digits = 4, width = 11),
        collapse = " "
      ))
    })
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
  tv <- do.call(rbind, lapply(seq_len(n), function(i) {
    unlist(fun(thetas, covs[i, , drop = FALSE]))
  }))

  tab <- covs
  etas <- stats::rnorm(n, sd = 0.2)
  for (p in out$functionListName) {
    idx <- if (p %in% names(out$etaMap)) out$etaMap[[p]] else NULL # secondaries have none
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
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(makeTable(out, thetas), f)

  v <- verifyParamFunction(out, f, thetas, quiet = TRUE)

  # scalar TRUE/FALSE, usable directly in an if
  expect_length(as.logical(v), 1L)
  expect_true(as.logical(v))
  expect_true(if (v) TRUE else FALSE)

  d <- attr(v, "checks")
  expect_s3_class(d, "data.frame")
  expect_named(d, c(
    "PARAMETER", "TABLECOLUMN", "N", "MAXABSDIFF",
    "MAXRELDIFF", "PASS"
  ))
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
  broken <- eval(parse(text = sub(
    "thetas\\[2\\]", "thetas[3]",
    paste(out$code, collapse = "\n")
  )))
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

  expect_warning(
    v <- verifyParamFunction(out, f, thetas, quiet = TRUE),
    "not columns of"
  )
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
  expect_warning(
    v <- verifyParamFunction(out, f, thetas, quiet = TRUE),
    "not a column of"
  )
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
  tab$ETA3 <- NULL # CL is tabled but cannot be reduced to typical

  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(tab, f)
  expect_warning(
    v <- verifyParamFunction(out, f, thetas, quiet = TRUE),
    "Cannot recover typical values"
  )
  expect_true(is.na(attr(v, "checks")$PASS))
  expect_false(as.logical(v))
})

test_that("rows are deduplicated on the covariate combination", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  thetas <- run7Thetas()
  tab <- makeTable(out, thetas, n = 5)
  tab <- tab[rep(seq_len(nrow(tab)), each = 4), ] # 4 records per subject

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

  v <- verifyParamFunction(out, f, thetas, quiet = TRUE)
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

test_that("secondary parameters are skipped by default", {
  out <- suppressWarnings(
    createParamFunction(modFile,
      parameters = c("CL", "V"), quiet = TRUE,
      secondary = list(AUC = "80 / CL")
    )
  )
  thetas <- run7Thetas()
  f <- withr::local_tempfile(fileext = ".tab")
  # makeTable() tables every functionListName entry, so the file gets an AUC
  # column too - but the default check still skips it (secondaries are excluded
  # from `parameters` unless named).
  writeNmTable(makeTable(out, thetas), f)

  v <- verifyParamFunction(out, f, thetas, quiet = TRUE)
  d <- attr(v, "checks")
  expect_setequal(d$PARAMETER, c("CL", "V")) # AUC not checked
  expect_true(as.logical(v))

  # ... but naming it explicitly forces the check. AUC is in the table (a raw
  # column) yet has no TVAUC column and no exponential-IIV pattern, so it
  # cannot be reduced to a typical value -> warn, PASS = NA.
  expect_warning(
    v2 <- verifyParamFunction(out, f, thetas, parameters = "AUC", quiet = TRUE),
    "Cannot recover typical values"
  )
  expect_true(is.na(attr(v2, "checks")$PASS))
  expect_false(as.logical(v2))
})

test_that("naming a parameter the function does not return warns rather than crashing", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE)
  )
  thetas <- run7Thetas()
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(makeTable(out, thetas), f)

  # "Cmax" is neither a $PK parameter nor a defined secondary.
  expect_warning(
    v <- verifyParamFunction(out, f, thetas, parameters = "Cmax", quiet = TRUE),
    "not returned by the generated function"
  )
  d <- attr(v, "checks")
  expect_equal(d$PARAMETER, "Cmax")
  expect_true(is.na(d$PASS))
  expect_equal(d$N, 0L)
  expect_false(as.logical(v))
})

test_that("bad input is rejected", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", quiet = TRUE)
  )
  expect_error(
    verifyParamFunction(list(), "nowhere.tab", 1),
    "returned by createParamFunction"
  )
  expect_error(
    verifyParamFunction(out, "does-not-exist.tab", 1),
    "Table file not found"
  )
})

## --- table formats other than "banner + whitespace" ------------------------

writeCsvTable <- function(df, path, banner = FALSE) {
  lines <- c(
    paste(names(df), collapse = ","),
    apply(df, 1, function(r) paste(formatC(as.numeric(r), format = "E", digits = 4), collapse = ","))
  )
  writeLines(if (banner) c("TABLE NO.  1", lines) else lines, path)
  path
}

test_that("a comma-separated table with no banner is read, not mis-parsed", {
  out <- suppressWarnings(createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE))
  thetas <- run7Thetas()
  f <- withr::local_tempfile(fileext = ".csv")
  writeCsvTable(makeTable(out, thetas), f)

  v <- verifyParamFunction(out, f, thetas, quiet = TRUE)

  expect_true(as.logical(v))
  expect_true(all(attr(v, "checks")$PASS))
})

test_that("a whitespace table with no banner is read", {
  out <- suppressWarnings(createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE))
  thetas <- run7Thetas()
  df <- makeTable(out, thetas)
  f <- withr::local_tempfile(fileext = ".tab")
  writeLines(c(
    paste(names(df), collapse = " "),
    apply(df, 1, function(r) paste(formatC(as.numeric(r), format = "E", digits = 4), collapse = " "))
  ), f)

  v <- verifyParamFunction(out, f, thetas, quiet = TRUE)

  expect_true(as.logical(v))
})

## --- rows the table cannot reconcile ---------------------------------------

test_that("a missVal covariate row is dropped when $PK has no guard for it", {
  ## run7's WT reference comes from the normalisation constant in (WT/75), so
  ## $PK has no IF(WT.EQ.-99) guard. A table row carrying WT = -99 therefore
  ## cannot be reconciled: NONMEM computed the tabled value from the real
  ## weight, which the table no longer shows. Comparing it measures the table.
  out <- suppressWarnings(createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE))
  expect_false(grepl("explicit missing-value handling", out$covRef$WT$source))

  thetas <- run7Thetas()
  df <- makeTable(out, thetas)
  ## Blank one row's WT after its CL was computed - exactly what a
  ## post-processed table looks like.
  df$WT[1] <- out$missVal
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(df, f)

  expect_warning(
    v <- verifyParamFunction(out, f, thetas, quiet = TRUE),
    "left out of the comparison"
  )
  ## and it names only the covariate that actually held missVal
  w <- tryCatch(verifyParamFunction(out, f, thetas, quiet = TRUE),
    warning = conditionMessage
  )
  expect_match(w, "WT")
  ## FOOD is a covariate of CL and did *not* hold missVal, so naming it would
  ## mean the warning had gone back to listing every unguarded covariate.
  ## (SEX would prove nothing here - pruning removes it from covRef entirely.)
  expect_false(grepl("FOOD", w))
  ## The remaining rows are consistent, so the check still passes.
  expect_true(as.logical(v))
  expect_true(all(attr(v, "checks")$N < nrow(df)))
})

## --- TV<P> is a naming convention, not a fact ------------------------------

test_that("a TV<P> column that is not P's typical value is not used as one", {
  ## run7's $PK ends:
  ##   TVMAT = THETA(6)     MAT = TVMAT * EXP(ETA(5))
  ##   TVD1  = THETA(7)     D1  = MAT * (1 - TVD1)
  ## so TVMAT really is MAT's typical value, while TVD1 is a dimensionless
  ## fraction that D1 is computed *from*. Comparing D1 against a TVD1 column
  ## compares a duration with a fraction and fails whatever the translation.
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("MAT", "D1"), quiet = TRUE)
  )
  expect_equal(unname(out$tvMap[["MAT"]]), "TVMAT")
  expect_false("D1" %in% names(out$tvMap))

  thetas <- run7Thetas()
  df <- makeTable(out, thetas)
  ## The model's own TVD1 - a fraction, nothing like D1.
  df$TVD1 <- thetas[7]
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(df, f)

  expect_warning(
    v <- verifyParamFunction(out, f, thetas, quiet = TRUE),
    "TVD1"
  )
  d <- attr(v, "checks")
  ## MAT verifies through its genuine typical-value column.
  expect_true(d$PASS[d$PARAMETER == "MAT"])
  ## D1 is reported as not checkable, not as a failure.
  expect_true(is.na(d$PASS[d$PARAMETER == "D1"]))
})

test_that("a comparison-only covariate at missVal does not cost a row", {
  ## run7 tests GENO1 only as IF(GENO1.EQ.1), and its reference is 0. Neither
  ## -99 nor 0 is the tested level, so both take the same branch and the row
  ## reconciles - dropping it would lose coverage for nothing. This is the
  ## shape of a real model where -99 is simply another level, not a missing
  ## marker.
  out <- suppressWarnings(createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE))
  expect_false(out$covRef$GENO1$inArithmetic)
  expect_equal(out$covRef$GENO1$testedLevels, 1)

  thetas <- run7Thetas()
  df <- makeTable(out, thetas)
  df$GENO1[1] <- out$missVal
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(df, f)

  expect_no_warning(v <- verifyParamFunction(out, f, thetas, quiet = TRUE))
  expect_true(as.logical(v))
  expect_equal(attr(v, "checks")$N[1], nrow(df))
})

test_that("a covariate in arithmetic at missVal is dropped, and flagged as a model problem", {
  ## WT reaches (WT/75)**THETA, so -99 there is not something the model can
  ## have meant - NONMEM would raise a negative base to a fractional power.
  ## The row cannot be reconciled, and the control stream is worth a look.
  out <- suppressWarnings(createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE))
  expect_true(out$covRef$WT$inArithmetic)

  thetas <- run7Thetas()
  df <- makeTable(out, thetas)
  df$WT[1] <- out$missVal
  f <- withr::local_tempfile(fileext = ".tab")
  writeNmTable(df, f)

  w <- tryCatch(verifyParamFunction(out, f, thetas, quiet = TRUE),
    warning = conditionMessage
  )
  expect_match(w, "WT")
  expect_match(w, "arithmetic")
  expect_false(grepl("GENO1", w))

  ## and the row really is left out, not merely complained about
  v <- suppressWarnings(verifyParamFunction(out, f, thetas, quiet = TRUE))
  expect_true(all(attr(v, "checks")$N < nrow(df)))
})
