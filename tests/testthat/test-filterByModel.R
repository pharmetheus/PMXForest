modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")

simData <- function() {
  read.csv(system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv",
    package = "PMXForest"
  ))
}

## A minimal control stream with a given $INPUT and $DATA record.
tempMod <- function(input, dataRec) {
  f <- withr::local_tempfile(fileext = ".mod", .local_envir = parent.frame())
  writeLines(c(
    "$PROBLEM t", input, dataRec, "$PK", "CL = THETA(1)",
    "$THETA 1"
  ), f)
  f
}

test_that("run7's IGNOREs reproduce the records NONMEM actually used", {
  d <- simData()
  used <- filterByModel(d, modFile, quiet = TRUE)

  # xptab7 was written by this model and has 33885 data rows
  expect_equal(nrow(used), 33885)
  expect_equal(length(unique(used$ID)), 754)
  expect_equal(length(unique(d$ID)), 964)

  # the three conditions, applied by hand
  expect_equal(used, subset(d, TYPE != 2 & BLQ != 1 & ID != 895),
    ignore_attr = TRUE
  )
})

test_that("filtering changes the covariate quantiles it feeds", {
  d <- simData()
  used <- filterByModel(d, modFile, quiet = TRUE)
  expect_false(identical(
    getCovStats(d, "CRCL", idVar = "ID")$CRCL,
    getCovStats(used, "CRCL", idVar = "ID")$CRCL
  ))
})

test_that("columns are matched by position, not by name", {
  # run7's $INPUT calls position 9 ODV and position 10 DV; the csv header calls
  # them DV and LNDV. A name-based filter would read the wrong column.
  d <- simData()
  expect_equal(names(d)[9:10], c("DV", "LNDV"))

  renamed <- filterByModel(d, modFile, useInputNames = TRUE, quiet = TRUE)
  expect_equal(names(renamed)[9:10], c("ODV", "DV"))
  # the values are unchanged, only the labels move
  expect_equal(renamed$ODV, d$DV[d$TYPE != 2 & d$BLQ != 1 & d$ID != 895])
})

test_that("useInputNames returns only the columns the model reads", {
  d <- simData()
  mod <- nmReadModel(modFile)
  n <- length(nmInputPositions(mod)$names)
  expect_equal(ncol(filterByModel(d, modFile, useInputNames = TRUE, quiet = TRUE)), n)
  # the default keeps the caller's own frame intact
  expect_equal(names(filterByModel(d, modFile, quiet = TRUE)), names(d))
})

test_that("ACCEPT keeps only matching records", {
  f <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=@ ACCEPT=(WT.GT.70)")
  d <- data.frame(ID = 1:4, DV = 1, WT = c(60, 75, 80, 65))
  expect_equal(filterByModel(d, f, quiet = TRUE)$WT, c(75, 80))
})

test_that("a comma inside a list means OR", {
  f <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT.EQ.60,WT.EQ.80)")
  d <- data.frame(ID = 1:4, DV = 1, WT = c(60, 75, 80, 65))
  expect_equal(filterByModel(d, f, quiet = TRUE)$WT, c(75, 65))
})

test_that("several IGNORE statements are combined with OR", {
  f <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE(WT.EQ.60) IGNORE(ID.EQ.4)")
  d <- data.frame(ID = 1:4, DV = 1, WT = c(60, 75, 80, 65))
  expect_equal(filterByModel(d, f, quiet = TRUE)$ID, c(2, 3))
})

test_that("a bare = means equality and the = after the keyword is optional", {
  d <- data.frame(ID = 1:4, DV = 1, WT = c(60, 75, 80, 65))
  f1 <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT=60)")
  f2 <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE(WT.EQN.60)")
  expect_equal(filterByModel(d, f1, quiet = TRUE)$ID, c(2, 3, 4))
  expect_equal(filterByModel(d, f2, quiet = TRUE)$ID, c(2, 3, 4))
})

test_that("the full comparison operator family is understood", {
  d <- data.frame(ID = 1:5, DV = 1, WT = c(50, 60, 70, 80, 90))
  ge <- tempMod("$INPUT ID DV WT", "$DATA d.csv ACCEPT=(WT.GE.70)")
  le <- tempMod("$INPUT ID DV WT", "$DATA d.csv ACCEPT=(WT.LE.60)")
  ne <- tempMod("$INPUT ID DV WT", "$DATA d.csv ACCEPT=(WT.NE.70)")
  expect_equal(filterByModel(d, ge, quiet = TRUE)$WT, c(70, 80, 90))
  expect_equal(filterByModel(d, le, quiet = TRUE)$WT, c(50, 60))
  expect_equal(filterByModel(d, ne, quiet = TRUE)$WT, c(50, 60, 80, 90))
})

test_that(".EQN. and .NEN. behave as .EQ. and .NE.", {
  # run7's own $DATA uses .EQN., so the variant the reference model depends on
  # must be covered, and its partner with it.
  d <- data.frame(ID = 1:4, DV = 1, WT = c(50, 60, 60, 70))
  eqn <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE(WT.EQN.60)")
  nen <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE(WT.NEN.60)")
  eq <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE(WT.EQ.60)")
  ne <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE(WT.NE.60)")

  expect_equal(filterByModel(d, eqn, quiet = TRUE)$WT, c(50, 70))
  expect_equal(filterByModel(d, nen, quiet = TRUE)$WT, c(60, 60))
  # the .xxN. forms must agree with their plain partners
  expect_equal(
    filterByModel(d, eqn, quiet = TRUE), filterByModel(d, eq, quiet = TRUE)
  )
  expect_equal(
    filterByModel(d, nen, quiet = TRUE), filterByModel(d, ne, quiet = TRUE)
  )
})

test_that("a record the condition cannot be evaluated on is not selected", {
  # A missing value makes the comparison NA, and NA row-indexing keeps an
  # all-NA row rather than dropping it, so the result is forced to FALSE. That
  # reads the same way for both keywords: the condition did not fire. For
  # IGNORE the record is therefore kept, for ACCEPT it is dropped - which is
  # the safe direction in each case, but it is a silent decision about which
  # subjects reach the plot, so it is pinned here.
  d <- data.frame(ID = 1:4, DV = 1, WT = c(50, NA, 80, 90))

  ign <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT.GT.70)")
  acc <- tempMod("$INPUT ID DV WT", "$DATA d.csv ACCEPT=(WT.GT.70)")

  expect_equal(filterByModel(d, ign, quiet = TRUE)$ID, c(1L, 2L))
  expect_equal(filterByModel(d, acc, quiet = TRUE)$ID, c(3L, 4L))
  # no row is both kept by IGNORE and kept by ACCEPT except through the NA rule
  expect_false(2L %in% filterByModel(d, acc, quiet = TRUE)$ID)
})

test_that("several ACCEPT statements are alternatives, as several IGNOREs are", {
  # NM-TRAN: "Multiple IGNORE options with different lists may be used", the
  # conditions being joined by an implied .OR., and ACCEPT is "identical to the
  # IGNORE list option, except that it specifies conditions for acceptance".
  # So two ACCEPT statements accept the union, not the intersection.
  d <- data.frame(ID = 1:6, DV = 1, AGE = c(10, 20, 30, 40, 50, 60), SEX = c(1, 1, 2, 2, 1, 2))
  two <- tempMod("$INPUT ID DV AGE SEX", "$DATA d.csv ACCEPT=(AGE.GT.25) ACCEPT=(SEX.EQ.1)")
  one <- tempMod("$INPUT ID DV AGE SEX", "$DATA d.csv ACCEPT=(AGE.GT.25,SEX.EQ.1)")

  expect_equal(filterByModel(d, two, quiet = TRUE)$ID, 1:6)
  # a single list with both conditions must give the same answer
  expect_equal(
    filterByModel(d, two, quiet = TRUE), filterByModel(d, one, quiet = TRUE)
  )
  # and it is a union, not an intersection - the intersection is ID 5 alone
  expect_gt(nrow(filterByModel(d, two, quiet = TRUE)), 1)
})

test_that("a text value is refused with a message that names the problem", {
  # NONMEM allows IGNORE=(GEN='M'); the condition grammar shared with $PK has
  # no string literal, so it must say so rather than blaming $INPUT.
  d <- data.frame(ID = 1:2, DV = 1, GEN = c("M", "F"), stringsAsFactors = FALSE)
  f <- tempMod("$INPUT ID DV GEN", "$DATA d.csv IGNORE=(GEN.EQ.'M')")
  expect_error(filterByModel(d, f, quiet = TRUE), "compares against a text value")
  expect_error(filterByModel(d, f, quiet = TRUE), "numeric comparisons only")
})

test_that("the string family compares as text and the N family numerically", {
  # NM-TRAN: "With =, ==, /=, .EQ. and .NE., the value in the data record and
  # the value in the list are compared as character strings. Otherwise, they
  # are converted to numeric" - which is the case with .NEN. and .EQN.
  expect_equal(nmConditionToR("TYPE.EQ.2", "m.mod"), 'as.character(TYPE) == "2"')
  expect_equal(nmConditionToR("TYPE.NE.2", "m.mod"), 'as.character(TYPE) != "2"')
  expect_equal(nmConditionToR("OCC=1", "m.mod"), 'as.character(OCC) == "1"')
  expect_equal(nmConditionToR("OCC/=1", "m.mod"), 'as.character(OCC) != "1"')
  # the N variants and every inequality stay numeric
  expect_equal(nmConditionToR("TYPE.EQN.2", "m.mod"), "TYPE == 2")
  expect_equal(nmConditionToR("TYPE.NEN.2", "m.mod"), "TYPE != 2")
  expect_equal(nmConditionToR("WT.GT.70", "m.mod"), "WT > 70")
  # a bare = must not be taken out of >= or <=
  expect_equal(nmConditionToR("WT.GE.70", "m.mod"), "WT >= 70")
  expect_equal(nmConditionToR("WT>=70", "m.mod"), "WT >= 70")
  expect_equal(nmConditionToR("WT<=70", "m.mod"), "WT <= 70")
})

test_that("a text comparison that finds nothing numerically would is warned about", {
  # The table-file case the NM-TRAN help calls out: an integer 1 in the data is
  # written 1.0000E+00 in a table file, so IGNORE=(OCC.EQ.1) matches nothing.
  # We cannot see the file's text from a data frame, so the disagreement
  # between the two readings is the only available signal - and it is loud.
  d <- data.frame(ID = 1:4, DV = 1, OCC = c(1, 2, 1, 2))
  trap <- tempMod("$INPUT ID DV OCC", "$DATA d.csv IGNORE=(OCC.EQ.1.0000E+00)")

  expect_warning(
    out <- filterByModel(d, trap, quiet = TRUE),
    "selects no record"
  )
  expect_warning(filterByModel(d, trap, quiet = TRUE), "EQN")
  # nothing was dropped, because the text did not match
  expect_equal(nrow(out), 4L)

  # the numeric operator does the job and says nothing
  ok <- tempMod("$INPUT ID DV OCC", "$DATA d.csv IGNORE=(OCC.EQN.1.0000E+00)")
  expect_no_warning(out2 <- filterByModel(d, ok, quiet = TRUE))
  expect_equal(out2$ID, c(2L, 4L))

  # and a text comparison that does match is silent
  plain <- tempMod("$INPUT ID DV OCC", "$DATA d.csv IGNORE=(OCC.EQ.1)")
  expect_no_warning(out3 <- filterByModel(d, plain, quiet = TRUE))
  expect_equal(out3$ID, c(2L, 4L))
})

test_that("DROP columns still occupy a position", {
  f <- tempMod("$INPUT ID JUNK=DROP WT", "$DATA d.csv IGNORE=(WT.LT.70)")
  d <- data.frame(ID = 1:3, JUNK = 9, WT = c(60, 75, 80))
  expect_equal(filterByModel(d, f, quiet = TRUE)$WT, c(75, 80))
})

test_that("a SYNONYM=REAL pair can be referred to by either name", {
  d <- data.frame(ID = 1:3, CONC = c(1, 2, 3), WT = 70)
  f1 <- tempMod("$INPUT ID CONC=DV WT", "$DATA d.csv IGNORE=(DV.EQ.2)")
  f2 <- tempMod("$INPUT ID CONC=DV WT", "$DATA d.csv IGNORE=(CONC.EQ.2)")
  expect_equal(filterByModel(d, f1, quiet = TRUE)$CONC, c(1, 3))
  expect_equal(filterByModel(d, f2, quiet = TRUE)$CONC, c(1, 3))
})

test_that("a model with no filter returns the data unchanged", {
  f <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=@")
  d <- data.frame(ID = 1:3, DV = 1, WT = c(60, 75, 80))
  expect_equal(filterByModel(d, f, quiet = TRUE), d)
})

test_that("the conventional header markers do not warn, other characters do", {
  d <- data.frame(ID = 1:3, DV = 1, WT = c(60, 75, 80))
  fAt <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=@")
  fHash <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=#")
  fC <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=C")
  expect_silent(filterByModel(d, fAt, quiet = TRUE))
  expect_silent(filterByModel(d, fHash, quiet = TRUE))
  expect_warning(filterByModel(d, fC, quiet = TRUE), "IGNORE=C")
})

test_that("ACCEPT and IGNORE lists together are refused, as in NONMEM", {
  f <- tempMod("$INPUT ID DV WT", "$DATA d.csv ACCEPT=(WT.GT.70) IGNORE=(ID.EQ.1)")
  d <- data.frame(ID = 1:3, DV = 1, WT = c(60, 75, 80))
  expect_error(filterByModel(d, f, quiet = TRUE), "cannot both appear")
})

test_that("structural problems are reported clearly", {
  d <- data.frame(ID = 1:3, DV = 1, WT = c(60, 75, 80))

  tooFew <- tempMod("$INPUT ID DV WT AGE SEX", "$DATA d.csv IGNORE=(WT.LT.70)")
  expect_error(filterByModel(d, tooFew, quiet = TRUE), "at least the columns")

  noInput <- withr::local_tempfile(fileext = ".mod")
  writeLines(c("$PROBLEM t", "$DATA d.csv", "$PK", "CL = THETA(1)"), noInput)
  expect_error(filterByModel(d, noInput, quiet = TRUE), "No \\$INPUT record")

  unknown <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(NOPE.EQ.1)")
  expect_error(filterByModel(d, unknown, quiet = TRUE), "\\$INPUT does not declare")
})

test_that("a '.' column works under both comparison families", {
  d <- data.frame(ID = 1:3, DV = 1, WT = c("60", ".", "80"), stringsAsFactors = FALSE)

  ## .EQ. is a text comparison in NONMEM, so the column is compared as written.
  f <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT.EQ.60)")
  expect_equal(filterByModel(d, f, quiet = TRUE)$WT, c(".", "80"))

  ## .EQN. and the ordering comparisons read "." as 0, as NM-TRAN does.
  g <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT.EQN.60)")
  expect_equal(nrow(filterByModel(d, g, quiet = TRUE)), 2L)

  h <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT.GT.60)")
  expect_equal(nrow(filterByModel(d, h, quiet = TRUE)), 2L)

  ## and it says so, since this is a reading of the data rather than a literal
  expect_message(filterByModel(d, h), "placeholder")

  ## The value has to be 0 specifically. Every comparison above is false at 0
  ## AND at NA - an unresolvable row is simply not selected - so none of them
  ## can tell the two apart. A comparison that is TRUE at zero can.
  k <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT.LT.1)")
  kept <- filterByModel(d, k, quiet = TRUE)
  expect_equal(nrow(kept), 2L)
  expect_false("." %in% kept$WT)
})

test_that("extra columns beyond $INPUT are reported but harmless", {
  f <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT.LT.70)")
  d <- data.frame(ID = 1:3, DV = 1, WT = c(60, 75, 80), EXTRA = "x")
  expect_message(filterByModel(d, f), "beyond that are not read")
  expect_equal(filterByModel(d, f, quiet = TRUE)$WT, c(75, 80))
})

test_that("quiet = FALSE reports the condition and what it removed", {
  f <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT.LT.70)")
  d <- data.frame(ID = 1:3, DV = 1, WT = c(60, 75, 80))
  expect_message(filterByModel(d, f), "Applying IGNORE")
  expect_message(filterByModel(d, f), "Removed 1 of 3 record")
})

## --- against NONMEM's own account of what it used ------------------------
##
## The .lst prints how many records, subjects and observations NONMEM read
## after applying $DATA. That is an independent oracle - it comes from NONMEM,
## not from us - and it is the only check here that can catch filterByModel()
## agreeing with a hand-written subset while both are wrong.

## What NONMEM printed about the data it used.
lstCounts <- function(lstFile) {
  L <- readLines(lstFile, warn = FALSE)
  num <- function(pat) {
    hit <- grep(pat, L, value = TRUE)
    if (!length(hit)) {
      return(NA_integer_)
    }
    suppressWarnings(as.integer(sub(".*?:\\s*([0-9]+).*", "\\1", hit[1])))
  }
  list(
    records = num("NO\\. OF DATA RECS IN DATA SET:"),
    obs = num("TOT\\. NO\\. OF OBS RECS:"),
    ids = num("TOT\\. NO\\. OF INDIVIDUALS:")
  )
}

## NONMEM counts an observation as MDV == 0, or EVID == 0 where the data set
## carries no MDV column.
obsCount <- function(d) {
  if ("MDV" %in% names(d)) sum(d$MDV == 0) else sum(d$EVID == 0)
}

test_that("filterByModel reproduces the counts NONMEM printed in the .lst", {
  cases <- list(
    list(
      name = "run7", mod = "SimVal/run7.mod", lst = "SimVal/run7.lst",
      data = "SimVal/DAT-1-MI-PMX-2.csv", sep = ","
    ),
    ## A different shape entirely: ADVAN13, no MDV column, so the observation
    ## count goes through the EVID fallback.
    list(
      name = "tte_weibull", mod = "tte/tte_weibull.mod",
      lst = "tte/tte_weibull.lst", data = "tte/tte_data1.dat", sep = ""
    )
  )

  for (cs in cases) {
    mod <- system.file("extdata", cs$mod, package = "PMXForest")
    lst <- system.file("extdata", cs$lst, package = "PMXForest")
    dat <- system.file("extdata", cs$data, package = "PMXForest")
    skip_if(!all(nzchar(c(mod, lst, dat))), paste(cs$name, "not bundled"))

    want <- lstCounts(lst)
    expect_false(is.na(want$records), info = cs$name)

    raw <- utils::read.table(dat, sep = cs$sep, header = TRUE, comment.char = "")
    used <- filterByModel(raw, mod, quiet = TRUE)

    expect_equal(nrow(used), want$records, info = cs$name)
    expect_equal(length(unique(used$ID)), want$ids, info = cs$name)
    expect_equal(obsCount(used), want$obs, info = cs$name)
  }
})

test_that("the run7 case is an informative test of the filter", {
  ## A model whose $DATA removes nothing would pass the check above while
  ## proving nothing about filterByModel(). Assert that this one does remove
  ## something, so the test cannot quietly become vacuous if the bundled data
  ## is ever replaced.
  lst <- system.file("extdata", "SimVal/run7.lst", package = "PMXForest")
  raw <- simData()
  expect_gt(nrow(raw), lstCounts(lst)$records)
})

test_that("a text column is fine where the comparison is textual", {
  ## A comment column written as "." placeholders with the odd digit is
  ## ordinary, and NONMEM compares .EQ. as text, so it filters perfectly well.
  ## Rejecting the column outright contradicted that support and refused the
  ## model.
  d <- data.frame(
    C = c(".", ".", "7", ".", "2"),
    ID = 1:5, DV = c(1, 2, 3, 4, 5),
    stringsAsFactors = FALSE
  )
  mod <- tempMod("$INPUT C ID DV", "$DATA d.csv IGNORE=@ IGNORE(C.EQ.2)")
  expect_true(is.character(d$C))

  kept <- filterByModel(d, mod, quiet = TRUE)
  ## NONMEM compares the text "2", so only that row goes
  expect_equal(nrow(kept), 4L)
  expect_false("2" %in% kept$C)
  expect_true("7" %in% kept$C)
})

test_that("a '.' column answers a numeric comparison as 0", {
  d <- data.frame(
    C = c(".", ".", "7", ".", "2"),
    ID = 1:5, DV = c(1, 2, 3, 4, 5),
    stringsAsFactors = FALSE
  )
  mod <- tempMod("$INPUT C ID DV", "$DATA d.csv IGNORE=@ IGNORE(C.EQN.2)")
  kept <- filterByModel(d, mod, quiet = TRUE)
  expect_equal(nrow(kept), 4L)
  expect_false("2" %in% kept$C)
})

test_that("only the columns used numerically are complained about", {
  d <- data.frame(
    C = c(".", "2", "."), BLQ = c(".", ".", "1"),
    ID = 1:3, DV = c(1, 2, 3), stringsAsFactors = FALSE
  )
  mod <- tempMod(
    "$INPUT C BLQ ID DV",
    "$DATA d.csv IGNORE=@ IGNORE(C.EQ.2) IGNORE(BLQ.EQN.1)"
  )
  ## both columns are "." placeholders, so both resolve - C as text, BLQ as 0
  kept <- filterByModel(d, mod, quiet = TRUE)
  expect_equal(nrow(kept), 1L)
})

test_that("a '.' placeholder is read as 0 in a numeric comparison", {
  ## NONMEM data files carry "." wherever a field does not apply, and NM-TRAN
  ## reads it as 0 - a data set full of them is ordinary, not broken. Refusing
  ## the column meant refusing models NONMEM runs.
  d <- data.frame(
    DVID = c("1", ".", "7", "2", "."),
    ID = 1:5, DV = 1, stringsAsFactors = FALSE
  )
  mod <- tempMod("$INPUT DVID ID DV", "$DATA d.csv IGNORE=@ IGNORE=(DVID.GT.6)")
  kept <- filterByModel(d, mod, quiet = TRUE)
  ## only DVID = 7 exceeds 6; the "." rows are 0 and stay
  expect_equal(nrow(kept), 4L)
  expect_false("7" %in% kept$DVID)
  expect_equal(sum(kept$DVID == "."), 2L)

  ## and 0 rather than NA: IGNORE=(DVID.EQN.0) selects the placeholders only if
  ## they really became zero. With NA they would not be selected at all, and
  ## the row count above would be identical - which is why it proves nothing
  ## on its own.
  zero <- tempMod("$INPUT DVID ID DV", "$DATA d.csv IGNORE=@ IGNORE=(DVID.EQN.0)")
  expect_equal(nrow(filterByModel(d, zero, quiet = TRUE)), 3L)
})

test_that("real text in a numeric comparison is still refused", {
  ## "." is a convention. "unknown" is a data problem, and coercing it would
  ## drop rows the model kept without saying so.
  d <- data.frame(
    DVID = c("1", "unknown", "7"), ID = 1:3, DV = 1,
    stringsAsFactors = FALSE
  )
  mod <- tempMod("$INPUT DVID ID DV", "$DATA d.csv IGNORE=@ IGNORE=(DVID.GT.6)")
  expect_error(filterByModel(d, mod, quiet = TRUE), "read as text")
})

test_that("a column used only textually is left alone", {
  ## No coercion where none is needed - the text comparison wants the string.
  d <- data.frame(
    C = c(".", "2", "."), ID = 1:3, DV = 1,
    stringsAsFactors = FALSE
  )
  mod <- tempMod("$INPUT C ID DV", "$DATA d.csv IGNORE=@ IGNORE(C.EQ.2)")
  kept <- filterByModel(d, mod, quiet = TRUE)
  expect_equal(kept$C, c(".", "."))
})

test_that("coercion for a numeric comparison does not change a textual one", {
  ## A column used by both families must be read both ways: NONMEM compares
  ## .EQ. against the characters in the file, so "1.0" does not equal 1 and
  ## the record is kept. Coercing the whole column for the sake of .GT. made
  ## as.character(as.numeric("1.0")) == "1" and dropped it - and silenced the
  ## warning that exists to catch exactly that confusion.
  d <- data.frame(
    FLAG = c("1.0", "7.0", ".", "1"), ID = 1:4, DV = 1,
    stringsAsFactors = FALSE
  )
  textOnly <- tempMod("$INPUT FLAG ID DV", "$DATA d.csv IGNORE=(FLAG.EQ.1)")
  both <- tempMod("$INPUT FLAG ID DV", "$DATA d.csv IGNORE=(FLAG.EQ.1,FLAG.GT.5)")

  keptText <- filterByModel(d, textOnly, quiet = TRUE)
  keptBoth <- suppressWarnings(filterByModel(d, both, quiet = TRUE))

  ## "1.0" survives the .EQ. either way; only "1" and "7.0" go
  expect_true("1.0" %in% keptText$FLAG)
  expect_true("1.0" %in% keptBoth$FLAG)
  expect_false("1" %in% keptBoth$FLAG)
  expect_false("7.0" %in% keptBoth$FLAG)
})

## Subjects NONMEM reads but that carry no observation record ------------------

## run7's dose-only subjects, found by hand: after the three IGNOREs, every
## record they have left is a dose (EVID 1 or 4).
noObsIds <- function() {
  u <- subset(simData(), TYPE != 2 & BLQ != 1 & ID != 895)
  hasObs <- tapply(u$EVID == 0, u$ID, any)
  as.numeric(names(hasObs)[!hasObs])
}

## A minimal $PRED control stream.
predMod <- function(input) {
  f <- withr::local_tempfile(fileext = ".mod", .local_envir = parent.frame())
  writeLines(c(
    "$PROBLEM t", input, "$DATA d.csv IGNORE=@", "$PRED",
    "Y = THETA(1) + EPS(1)", "$THETA 1"
  ), f)
  f
}

test_that("dropNoObs removes run7's subjects that have only dose records", {
  d <- simData()
  ids <- noObsIds()
  expect_length(ids, 32)

  used <- filterByModel(d, modFile, dropNoObs = TRUE, quiet = TRUE)
  expect_equal(length(unique(used$ID)), 754 - 32)
  expect_equal(nrow(used), 33885 - 585)
  expect_equal(used, subset(d, TYPE != 2 & BLQ != 1 & ID != 895 & !ID %in% ids),
    ignore_attr = TRUE
  )
})

test_that("dropNoObs is off by default, so the subjects match NONMEM's count", {
  used <- filterByModel(simData(), modFile, quiet = TRUE)
  expect_equal(length(unique(used$ID)), 754)
  expect_true(all(noObsIds() %in% used$ID))
})

test_that("the report says how many subjects had no observation and how they were found", {
  expect_message(
    filterByModel(simData(), modFile, dropNoObs = TRUE),
    "32 subject\\(s\\) with no observation record.*EVID"
  )
})

test_that("without EVID or MDV, doses are recognised from AMT and RATE", {
  ## run7 with EVID dropped: NONMEM never sees it and NM-TRAN derives it from
  ## the dose items. In run7 a record is a dose exactly when AMT or RATE is
  ## non-zero, so the same 32 subjects must be found.
  src <- readLines(modFile)
  i <- grep("^\\$INPUT", src)
  expect_match(src[i], " EVID ")
  src[i] <- sub(" EVID ", " EVID=DROP ", src[i])
  f <- withr::local_tempfile(fileext = ".mod")
  writeLines(src, f)

  d <- simData()
  used <- filterByModel(d, f, dropNoObs = TRUE, quiet = TRUE)
  kept <- unique(subset(d, TYPE != 2 & BLQ != 1 & ID != 895)$ID)
  expect_equal(sort(setdiff(kept, used$ID)), sort(noObsIds()))
})

test_that("a steady-state infusion with AMT = 0 is still a dose", {
  ## Subject 2's only record is an SS infusion: AMT 0, RATE and SS non-zero.
  ## Looking at AMT alone would count it as an observation.
  f <- tempMod("$INPUT ID TIME AMT RATE SS DV", "$DATA d.csv IGNORE=@")
  d <- data.frame(
    ID = c(1, 1, 2), TIME = c(0, 1, 0), AMT = 0,
    RATE = c(10, 0, 10), SS = c(1, 0, 1), DV = c(0, 5, 0)
  )
  expect_equal(unique(filterByModel(d, f, dropNoObs = TRUE, quiet = TRUE)$ID), 1)
})

test_that("MDV decides where present, and MDV = 100 is not an observation", {
  ## Subject 2 has an EVID = 0 record with MDV = 1: an observation event with
  ## no observation. Subject 3's only candidate has MDV = 100, which NONMEM
  ## ignores during estimation.
  f <- tempMod("$INPUT ID TIME AMT EVID MDV DV", "$DATA d.csv IGNORE=@")
  d <- data.frame(
    ID = c(1, 1, 2, 2, 3, 3), TIME = c(0, 1, 0, 1, 0, 1),
    AMT = c(100, 0, 100, 0, 100, 0), EVID = c(1, 0, 1, 0, 1, 0),
    MDV = c(1, 0, 1, 1, 1, 100), DV = 0
  )
  expect_equal(unique(filterByModel(d, f, dropNoObs = TRUE, quiet = TRUE)$ID), 1)
})

test_that("a column dropped from $INPUT is not used", {
  ## NONMEM never sees a =DROP column, so MDV here says nothing: NM-TRAN
  ## derives it from AMT, and subject 1's second record is an observation.
  f <- tempMod("$INPUT ID TIME AMT MDV=DROP DV", "$DATA d.csv IGNORE=@")
  d <- data.frame(
    ID = c(1, 1, 2), TIME = c(0, 1, 0), AMT = c(100, 0, 100),
    MDV = 1, DV = 0
  )
  expect_equal(unique(filterByModel(d, f, dropNoObs = TRUE, quiet = TRUE)$ID), 1)
})

test_that("a $PRED model has only MDV to go by", {
  ## AMT and EVID are PREDPP items; in a $PRED model they are ordinary columns.
  d <- data.frame(
    ID = c(1, 1, 2, 2), AMT = c(100, 0, 100, 100),
    EVID = c(1, 0, 1, 1), DV = 1
  )
  ## No MDV: every record is an observation, so nobody goes
  f <- predMod("$INPUT ID AMT EVID DV")
  expect_equal(nrow(filterByModel(d, f, dropNoObs = TRUE, quiet = TRUE)), 4)

  d$MDV <- c(0, 0, 1, 1)
  f2 <- predMod("$INPUT ID AMT EVID DV MDV")
  expect_equal(unique(filterByModel(d, f2, dropNoObs = TRUE, quiet = TRUE)$ID), 1)
})

test_that("'.' in a dose column reads as 0, as NM-TRAN reads it", {
  f <- tempMod("$INPUT ID TIME AMT DV", "$DATA d.csv IGNORE=@")
  d <- data.frame(
    ID = c(1, 1, 2), TIME = c(0, 1, 0), AMT = c("100", ".", "100"), DV = 0,
    stringsAsFactors = FALSE
  )
  expect_equal(unique(filterByModel(d, f, dropNoObs = TRUE, quiet = TRUE)$ID), 1)
})

test_that("a subject is a contiguous block of ID, as NONMEM reads it", {
  ## ID 1 comes back after ID 2. NONMEM starts a new individual there, and
  ## that one holds only a dose.
  f <- tempMod("$INPUT ID TIME AMT DV", "$DATA d.csv IGNORE=@")
  d <- data.frame(
    ID = c(1, 1, 2, 1), TIME = c(0, 1, 0, 5), AMT = c(100, 0, 0, 100), DV = 0
  )
  used <- filterByModel(d, f, dropNoObs = TRUE, quiet = TRUE)
  expect_equal(used$TIME, c(0, 1, 0))
})

test_that("dropNoObs needs an ID item", {
  f <- tempMod("$INPUT TIME AMT DV", "$DATA d.csv IGNORE=@")
  expect_error(
    filterByModel(data.frame(TIME = 0, AMT = 0, DV = 1), f,
      dropNoObs = TRUE, quiet = TRUE
    ),
    "ID"
  )
})
