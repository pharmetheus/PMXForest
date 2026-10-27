modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")
lstFile <- system.file("extdata", "SimVal/run7.lst", package = "PMXForest")
datFile <- system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")

test_that("it reproduces NONMEM's own counts and says so", {
  v <- verifyFilterByModel(modFile, quiet = TRUE)

  expect_length(as.logical(v), 1L)
  expect_true(as.logical(v))
  expect_true(if (v) TRUE else FALSE)

  d <- attr(v, "checks")
  expect_s3_class(d, "data.frame")
  expect_named(d, c("CHECK", "FILTERED", "NONMEM", "DIFF", "PASS"))
  expect_equal(d$CHECK, c("RECORDS", "SUBJECTS", "OBSERVATIONS"))
  expect_true(all(d$PASS))
  expect_equal(d$NONMEM, c(33885, 754, 6284))
  expect_true(all(d$DIFF == 0))
})

test_that("a data frame can be supplied instead of letting it find $DATA", {
  raw <- read.csv(datFile)
  v <- verifyFilterByModel(modFile, data = raw, lstFile = lstFile, quiet = TRUE)
  expect_true(as.logical(v))
})

test_that("it fails, rather than errors, when the data is not what was used", {
  raw <- read.csv(datFile)
  v <- suppressWarnings(
    verifyFilterByModel(modFile, data = raw[raw$ID > 100, ], lstFile = lstFile, quiet = TRUE)
  )
  expect_false(as.logical(v))
  expect_false(all(attr(v, "checks")$PASS))
})

test_that("it reports whether the test was worth anything", {
  ## A $DATA that removes nothing agrees with the .lst trivially.
  v <- verifyFilterByModel(modFile, quiet = TRUE)
  expect_true(attr(v, "informative"))
  expect_gt(attr(v, "removed"), 0)
})

test_that("a model whose $DATA filters nothing is flagged as uninformative", {
  mod <- system.file("extdata", "tte/tte_weibull.mod", package = "PMXForest")
  skip_if(!nzchar(mod), "tte model not bundled")
  v <- verifyFilterByModel(mod, quiet = TRUE)
  expect_true(as.logical(v))
  expect_false(attr(v, "informative"))
})

test_that("it says what it needs when it cannot find the pieces", {
  f <- withr::local_tempfile(fileext = ".mod")
  writeLines(c(
    "$PROBLEM t", "$INPUT ID TIME DV", "$DATA nowhere.csv IGNORE=@",
    "$PK", "CL = THETA(1)", "$THETA 1"
  ), f)
  expect_error(verifyFilterByModel(f, quiet = TRUE), "\\.lst")
})

test_that("it warns when the .lst describes a different control stream", {
  ## The counts only describe the model NONMEM actually read. If $INPUT has
  ## moved since, the columns no longer line up and a comparison is noise.
  f <- withr::local_tempfile(fileext = ".mod")
  writeLines(sub("^\\$INPUT      NO ID", "$INPUT      ID NO", readLines(modFile)), f)
  file.copy(lstFile, sub("\\.mod$", ".lst", f))
  expect_warning(
    verifyFilterByModel(f, data = read.csv(datFile), quiet = TRUE),
    "INPUT"
  )
})

test_that("a $INPUT differing only in DROP items is caught as stale", {
  ## A .mod edited after the run to add (or remove) a =DROP column shifts every
  ## position after it, because a DROP item still occupies a data column - the
  ## NONMEM guide's own example has DAT1=DROP as the second column of the file.
  ## Comparing names alone misses this entirely: DROP items are stripped from
  ## the name list, so both sides look the same while the mapping has moved.
  raw <- read.csv(datFile)

  ## a control stream carrying three extra DROP items the run did not have
  edited <- withr::local_tempfile(fileext = ".mod")
  L <- readLines(modFile)
  i <- grep("^\\$DATA", L)[1]
  writeLines(append(L, "            ADTM=DROP DATETIME=DROP USUBJID=DROP", after = i - 1L), edited)
  file.copy(lstFile, sub("\\.mod$", ".lst", edited))

  ## The warning comes first; filterByModel() then refuses on the column count,
  ## which is the same order a real stale model produces. Both are wanted - the
  ## warning is what explains the refusal.
  expect_warning(
    try(verifyFilterByModel(edited, data = raw, quiet = TRUE), silent = TRUE),
    "INPUT"
  )
  w <- tryCatch(
    withCallingHandlers(
      try(verifyFilterByModel(edited, data = raw, quiet = TRUE), silent = TRUE),
      warning = function(x) stop(conditionMessage(x))
    ),
    error = conditionMessage
  )
  expect_match(w, "declares 42 column\\(s\\) where the run read 39")
  expect_match(w, "=DROP item still occupies a column")
})

test_that("an unedited model raises no staleness warning", {
  raw <- read.csv(datFile)
  expect_no_warning(verifyFilterByModel(modFile, data = raw, quiet = TRUE))
})

test_that("the echoed control stream is found without an NM-TRAN MESSAGES marker", {
  ## Not every .lst carries that line - one real run does not - and the old
  ## code fell back to a fixed 400-line window from $PROBLEM. A control stream
  ## longer than that would have been truncated mid-echo and the $INPUT
  ## comparison would have run on half a record, silently.
  raw <- read.csv(datFile)
  edited <- withr::local_tempfile(fileext = ".mod")
  L <- readLines(modFile)
  i <- grep("^\\$DATA", L)[1]
  writeLines(append(L, "            ADTM=DROP DATETIME=DROP", after = i - 1L), edited)

  ## a .lst whose echo is followed by output, with no NM-TRAN MESSAGES line,
  ## and padded past 400 lines so the old fallback would have missed the end
  lstL <- readLines(lstFile)
  lstL <- lstL[!grepl("NM-TRAN MESSAGES", lstL)]
  ## Pad the echo itself past 400 lines. Without this the old fixed window
  ## still covered the whole control stream, so the test would pass under the
  ## very code it was written to condemn.
  at <- grep("^\\s*\\$EST", lstL)[1]
  lstL <- append(lstL, rep(";; padding", 400L), after = at - 1L)
  writeLines(lstL, sub("\\.mod$", ".lst", edited))

  expect_warning(
    try(verifyFilterByModel(edited, data = raw, quiet = TRUE), silent = TRUE),
    "INPUT"
  )
})

test_that("a .lst with no recognisable echo is not treated as agreement", {
  ## Silence has to mean "checked and matched", never "could not look".
  raw <- read.csv(datFile)
  edited <- withr::local_tempfile(fileext = ".mod")
  file.copy(modFile, edited, overwrite = TRUE)
  lstL <- readLines(lstFile)
  lstL <- lstL[!grepl("^\\s*\\$", lstL)]
  writeLines(lstL, sub("\\.mod$", ".lst", edited))

  expect_warning(
    verifyFilterByModel(edited, data = raw, quiet = TRUE),
    "could not be read"
  )
})

test_that("a data file short of the $INPUT the run read is named as such", {
  ## The useful diagnosis is not "your data frame is too narrow". It is that
  ## the model and the run agree on N columns and this file has fewer, so it
  ## is not the data the run used - every position after the gap is shifted.
  raw <- read.csv(datFile)
  short <- raw[, seq_len(ncol(raw) - 3L)]

  w <- tryCatch(
    verifyFilterByModel(modFile, data = short, lstFile = lstFile, quiet = TRUE),
    error = conditionMessage
  )
  expect_match(w, "39")
  expect_match(w, as.character(ncol(short)))
  expect_match(w, "not the data")
})

test_that("rawRows counts the data as read, not as reconstructed from a pass", {
  raw <- read.csv(datFile)

  v <- verifyFilterByModel(modFile, data = raw, lstFile = lstFile, quiet = TRUE)
  expect_equal(attr(v, "rawRows"), nrow(raw))
  ## `removed` is measured against NONMEM's own count, so it is independent of
  ## what filterByModel() did - which is the whole point of the comparison.
  expect_equal(attr(v, "removed"), nrow(raw) - attr(v, "checks")$NONMEM[1])

  ## On a FAIL the two no longer reconcile, and that is exactly when a report
  ## must not quietly reconstruct the row count from the numbers that disagree.
  short <- raw[raw$ID > 100, ]
  f <- suppressWarnings(
    verifyFilterByModel(modFile, data = short, lstFile = lstFile, quiet = TRUE)
  )
  expect_false(as.logical(f))
  expect_equal(attr(f, "rawRows"), nrow(short))
})

test_that("the observation basis falls back from MDV to EVID to the dose items", {
  ## These branches decide what OBSERVATIONS is compared against, and a wrong
  ## basis makes the comparison silently meaningless rather than fail.
  obsOf <- function(d, input, pred = FALSE) {
    f <- withr::local_tempfile(fileext = ".mod")
    writeLines(c(
      "$PROBLEM t", input, "$DATA d.csv IGNORE=@",
      if (pred) "$PRED" else "$PK", "Y = THETA(1)", "$THETA 1"
    ), f)
    mod <- nmReadModel(f)
    r <- nmObsRecords(mod, d, nmInputPositions(mod))
    list(n = sum(r$obs), basis = r$basis)
  }
  d <- data.frame(
    ID = 1:4, MDV = c(0, 1, 0, 0), EVID = c(0, 1, 0, 2),
    AMT = c(0, 10, 0, 0), RATE = 0
  )
  expect_equal(obsOf(d, "$INPUT ID MDV EVID AMT RATE"), list(n = 3L, basis = "MDV"))
  expect_equal(obsOf(d, "$INPUT ID MDV=DROP EVID AMT RATE"), list(n = 2L, basis = "EVID"))
  expect_equal(
    obsOf(d, "$INPUT ID MDV=DROP EVID=DROP AMT RATE"),
    list(n = 3L, basis = "AMT/RATE")
  )
  expect_equal(
    obsOf(d, "$INPUT ID MDV=DROP EVID=DROP AMT=DROP RATE=DROP"),
    list(n = 4L, basis = "NONE")
  )
  ## a $PRED model has no EVID or dose items to fall back on
  expect_equal(
    obsOf(d, "$INPUT ID MDV=DROP EVID AMT RATE", pred = TRUE),
    list(n = 4L, basis = "NONE")
  )

  ## An empty field is a null data item, like ".", and NM-TRAN reads both as
  ## 0 ($DATA help, NULL=): read.csv() gives NA for ",,", so NA is MDV = 0.
  na <- data.frame(ID = 1:3, MDV = c(0, NA, 1))
  expect_equal(obsOf(na, "$INPUT ID MDV")$n, 2L)
})

test_that("the .lst is looked for under each name PsN and NONMEM use", {
  dir <- withr::local_tempdir()
  m <- file.path(dir, "run1.mod")
  file.create(m)

  expect_null(nmFindLst(m))

  ## .res, .out and PsN's NM_run1/psn.lst are all real layouts; a fallback that
  ## silently never fires would send the caller to "no .lst found" instead.
  for (nm in c("run1.res", "run1.out")) {
    f <- file.path(dir, nm)
    file.create(f)
    expect_equal(normalizePath(nmFindLst(m)), normalizePath(f))
    unlink(f)
  }
  dir.create(file.path(dir, "NM_run1"))
  psn <- file.path(dir, "NM_run1", "psn.lst")
  file.create(psn)
  expect_equal(normalizePath(nmFindLst(m)), normalizePath(psn))

  ## .lst beside the model wins over the PsN copy
  lst <- file.path(dir, "run1.lst")
  file.create(lst)
  expect_equal(normalizePath(nmFindLst(m)), normalizePath(lst))
})

test_that("a changed $DATA filter is reported even when $INPUT still matches", {
  ## The two drift checks are independent: editing only the filter leaves
  ## every column in place, so the $INPUT comparison passes and this is the
  ## only thing standing between the caller and counts for a different subset.
  dir <- withr::local_tempdir()
  edited <- file.path(dir, "run7.mod")
  L <- readLines(modFile)
  i <- grep("IGNORE\\(TYPE\\.EQN\\.2\\)", L)
  ## Not skip_if: a condition that stopped matching would turn this into a
  ## test that reports success while asserting nothing.
  expect_length(i, 1L)
  L[i] <- sub("IGNORE\\(TYPE\\.EQN\\.2\\)", "IGNORE(TYPE.EQN.3)", L[i])
  writeLines(L, edited)
  expect_true(any(grepl("TYPE\\.EQN\\.3", readLines(edited))))
  file.copy(lstFile, file.path(dir, "run7.lst"))

  expect_warning(
    nmCheckLstMatches(edited, file.path(dir, "run7.lst")),
    "\\$DATA filter"
  )
  ## and the unedited model does not warn
  writeLines(readLines(modFile), edited)
  expect_silent(nmCheckLstMatches(edited, file.path(dir, "run7.lst")))
})

test_that("print() reports the verdict, the count and an uninformative check", {
  v <- verifyFilterByModel(modFile, quiet = TRUE)
  out <- paste(capture.output(print(v)), collapse = "\n")
  expect_match(out, "PASS")
  expect_match(out, "3/3")
  expect_match(out, "record\\(s\\) removed")
  expect_match(out, "OBSERVATIONS")

  ## an uninformative result must say so instead of showing a count that
  ## would read as evidence
  attr(v, "informative") <- FALSE
  out2 <- paste(capture.output(print(v)), collapse = "\n")
  expect_match(out2, "proves nothing")
  expect_false(grepl("record\\(s\\) removed", out2))

  expect_identical(withVisible(print(v))$visible, FALSE)
})

test_that("observations are counted under the $INPUT names, not the file's", {
  ## Rename the file's EVID column. NONMEM reads it by position and so must
  ## the observation count; a count by the file's names finds no EVID at all.
  d <- read.csv(datFile)
  expect_equal(names(d)[11], "EVID")
  names(d)[11] <- "EVENT"
  v <- verifyFilterByModel(modFile, data = d, quiet = TRUE)
  chk <- attr(v, "checks")
  expect_equal(chk$FILTERED[chk$CHECK == "OBSERVATIONS"], 6284)
  expect_equal(attr(v, "obsBasis"), "EVID")
})
