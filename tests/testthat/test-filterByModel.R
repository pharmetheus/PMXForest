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

test_that("a text column used by the filter is refused rather than compared", {
  f <- tempMod("$INPUT ID DV WT", "$DATA d.csv IGNORE=(WT.EQ.60)")
  d <- data.frame(ID = 1:3, DV = 1, WT = c("60", ".", "80"))
  expect_error(filterByModel(d, f, quiet = TRUE), "read as text rather than numbers")
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
