modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")
tteFile <- system.file("extdata", "tte/tte_weibull.mod", package = "PMXForest")

## Final estimates for run7. The last rows of a .ext hold standard errors and
## other special records, so select the estimates row explicitly.
run7Thetas <- function() {
  ext <- getExt(system.file("extdata", "SimVal/run7.ext", package = "PMXForest"))
  as.numeric(ext[ext$ITERATION == -1000000000, 2:15])
}

## Write a minimal control stream and return its path.
tempMod <- function(lines) {
  f <- withr::local_tempfile(fileext = ".mod", .local_envir = parent.frame())
  writeLines(lines, f)
  f
}

test_that("the return value has the documented shape", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE)
  )
  expect_named(out, c(
    "code", "functionListName", "primaryNames",
    "secondaryNames", "noBaseThetas", "covRef",
    "etaMap", "tvMap", "modFile", "missVal"
  ))
  expect_s3_class(out$code, "pmxParamFunction")
  expect_type(out$code, "character")
  expect_equal(out$functionListName, c("CL", "V"))
  expect_equal(out$primaryNames, c("CL", "V"))
  expect_equal(out$secondaryNames, character(0))
  expect_equal(out$noBaseThetas, 14)
})

test_that("the generated source parses and is a conforming parameter function", {
  out <- suppressWarnings(createParamFunction(modFile, quiet = TRUE))
  fun <- eval(parse(text = out$code))
  expect_type(fun, "closure")
  expect_equal(names(formals(fun)), c("thetas", "df", "..."))
})

test_that("covariates are those in $INPUT that $PK reads before assigning", {
  out <- suppressWarnings(createParamFunction(modFile, quiet = TRUE))
  expect_setequal(
    names(out$covRef),
    c("SEX", "GENO4", "FORM", "FOOD", "WT", "GENO1", "GENO3")
  )
  # GENO2 is assigned a constant before it is read, so it is a local
  expect_false("GENO2" %in% names(out$covRef))
  # TVCL is assigned before it is read and is never a covariate
  expect_false("TVCL" %in% names(out$covRef))
})

test_that("a covariate re-imputed by IF(X.EQ.missVal) is still a covariate", {
  # Rule 1 of nmCovRef() reads exactly this idiom, but the old classification
  # excluded any $INPUT name $PK assigned - so the documented rule could never
  # fire, and the generated function referred to an unbound WT.
  f <- tempMod(c(
    "$PROBLEM rule 1",
    "$INPUT ID TIME DV AMT WT",
    "$DATA data.csv IGNORE=@",
    "$PK",
    "IF(WT.EQ.-99) WT = 75",
    "TVCL = THETA(1)",
    "CL = TVCL*(WT/75)**0.75",
    "V = THETA(2)",
    "$THETA (0,5) (0,50)"
  ))
  out <- createParamFunction(f, parameters = c("CL", "V"), quiet = TRUE)

  expect_equal(names(out$covRef), "WT")
  expect_equal(out$covRef$WT$value, 75)
  expect_true(any(grepl('df[["WT"]]', out$code, fixed = TRUE)))

  fun <- eval(parse(text = paste(out$code, collapse = "\n")))
  expect_equal(
    fun(thetas = c(5, 50), df = data.frame(WT = -99))$CL,
    fun(thetas = c(5, 50), df = data.frame(WT = 75))$CL
  )
  expect_gt(
    fun(thetas = c(5, 50), df = data.frame(WT = 90))$CL,
    fun(thetas = c(5, 50), df = data.frame(WT = 75))$CL
  )
})

test_that("a symbol read before it is bound is refused at generation time", {
  base <- c(
    "$PROBLEM unbound", "$INPUT ID TIME DV AMT WT",
    "$DATA data.csv IGNORE=@", "$PK"
  )
  tail <- c("V = THETA(2)", "$THETA (0,7) (0,3)")

  # not in $INPUT and never assigned: a NONMEM reserved variable
  expect_error(
    createParamFunction(
      tempMod(c(
        base,
        "IF(NEWIND.NE.2) CNT = 0", "TVCL = THETA(1)",
        "CL = TVCL*(WT/75)*CNT", tail
      )),
      parameters = c("CL", "V"), quiet = TRUE
    ),
    "NEWIND"
  )
  expect_error(
    createParamFunction(
      tempMod(c(
        base,
        "IF(NEWIND.NE.2) CNT = 0", "TVCL = THETA(1)",
        "CL = TVCL*(WT/75)*CNT", tail
      )),
      parameters = c("CL", "V"), quiet = TRUE
    ),
    "supplied by NONMEM"
  )
  # read above its own assignment: NONMEM would carry a value over from the
  # previous data record, which a one-row parameter function cannot do
  expect_error(
    createParamFunction(
      tempMod(c(
        base,
        "CL = TVCL*(WT/75)", "TVCL = THETA(1)", tail
      )),
      parameters = c("CL", "V"), quiet = TRUE
    ),
    "reads TVCL.*before assigning it"
  )
  # but an exhaustive scm-style branch chain, as run7.mod writes, is accepted:
  # FRELWT is assigned only inside IFs with no ELSE, which no static analysis
  # can tell apart from a genuinely non-exhaustive branch
  scm <- createParamFunction(
    tempMod(c(
      base,
      "TVCL = THETA(1)",
      "IF(WT.GT.100) FRELWT = 1",
      "IF(WT.LE.100) FRELWT = 2",
      "CL = TVCL*FRELWT*(WT/75)", tail
    )),
    parameters = c("CL", "V"), quiet = TRUE
  )
  expect_equal(names(scm$covRef), "WT")
  expect_false("FRELWT" %in% names(scm$covRef))
})

test_that("covRef is validated rather than spliced in unchecked", {
  expect_error(
    suppressWarnings(createParamFunction(modFile,
      parameters = "CL",
      covRef = list(WGT = 70), quiet = TRUE
    )),
    "does not use: WGT"
  )
  for (bad in list(
    list(WT = "seventy"), list(WT = c(70, 80)),
    list(WT = NA), list(WT = Inf)
  )) {
    expect_error(
      suppressWarnings(createParamFunction(modFile,
        parameters = "CL",
        covRef = bad, quiet = TRUE
      )),
      "single finite number"
    )
  }
  expect_error(
    suppressWarnings(createParamFunction(modFile,
      parameters = "CL",
      covRef = list(70), quiet = TRUE
    )),
    "must be named"
  )
  # a valid override still works and is recorded as supplied
  ok <- suppressWarnings(createParamFunction(modFile,
    parameters = "CL",
    covRef = list(WT = 70), quiet = TRUE
  ))
  expect_equal(ok$covRef$WT$value, 70)
  expect_match(ok$covRef$WT$source, "supplied through covRef")
})

test_that("covariate references are taken from the control stream", {
  out <- suppressWarnings(createParamFunction(modFile, quiet = TRUE))

  # Rule 2a: the ";  Most common" branch
  expect_equal(out$covRef$SEX$value, 1)
  expect_match(out$covRef$SEX$source, "Most common")
  expect_true(out$covRef$SEX$confident)
  expect_equal(out$covRef$FOOD$value, 1)
  expect_equal(out$covRef$GENO4$value, 0)

  # Rule 3: the normalisation constant of (WT/75)
  expect_equal(out$covRef$WT$value, 75)
  expect_match(out$covRef$WT$source, "normalisation constant")

  # Rule 4: proposed, and flagged as such
  expect_equal(out$covRef$GENO1$value, 0)
  expect_false(out$covRef$GENO1$confident)
})

test_that("an inferred reference warns", {
  expect_warning(
    createParamFunction(modFile, quiet = TRUE),
    "inferred rather than read"
  )
})

test_that("covRef overrides a derived reference and is recorded as supplied", {
  out <- suppressWarnings(
    createParamFunction(modFile, covRef = list(WT = 70), quiet = TRUE)
  )
  expect_equal(out$covRef$WT$value, 70)
  expect_match(out$covRef$WT$source, "covRef")
  expect_true(any(grepl("else 70", out$code)))
})

test_that("the covariate preamble reads df[[cov]], so it cannot partial-match", {
  # `$` partial-matches on a data frame: data.frame(WTKG = 90)$WT is 90. A
  # data set carrying a longer name from the same family but not the covariate
  # itself must fall back to the reference, not silently use the other column.
  f <- tempMod(c(
    "$PROBLEM partial matching",
    "$INPUT ID TIME DV AMT WT WTKG",
    "$DATA data.csv IGNORE=@",
    "$PK",
    "TVCL = THETA(1)",
    "CL = TVCL*(WT/75)**0.75",
    "V = THETA(2)",
    "$THETA (0,7) (0,3)"
  ))
  out <- createParamFunction(f, parameters = c("CL", "V"), quiet = TRUE)

  expect_true(any(grepl('df[["WT"]]', out$code, fixed = TRUE)))
  expect_false(any(grepl("df$WT", out$code, fixed = TRUE)))

  fun <- eval(parse(text = paste(out$code, collapse = "\n")))
  atRef <- fun(thetas = c(7, 3), df = data.frame(FOO = 1))$CL
  expect_equal(fun(thetas = c(7, 3), df = data.frame(WTKG = 90))$CL, atRef)
  expect_equal(fun(thetas = c(7, 3), df = data.frame(WT = 75))$CL, atRef)
  # a real WT is still read
  expect_false(isTRUE(all.equal(
    fun(thetas = c(7, 3), df = data.frame(WT = 90))$CL, atRef
  )))
})

test_that("a covariate with no derivable reference is refused, not guessed", {
  # EXPO enters tte_weibull.mod as THETA(4)*EXPO: no branch, no normalisation.
  expect_error(
    createParamFunction(tteFile, quiet = TRUE),
    "No reference value could be derived"
  )
  expect_error(createParamFunction(tteFile, quiet = TRUE), "EXPO")
  # ... and supplying it is enough to proceed
  out <- createParamFunction(tteFile, covRef = list(EXPO = 0), quiet = TRUE)
  expect_equal(out$noBaseThetas, 4)
  expect_equal(out$covRef$AGE$value, 50)
})

test_that("generated values match the hand-written function from the Walkthrough", {
  handWritten <- function(thetas, df, ...) {
    CLFOOD <- 1
    if (any(names(df) == "FOOD") && df$FOOD != -99 && df$FOOD == 0) CLFOOD <- 1 + thetas[11]
    FRELFORM <- 1
    if (any(names(df) == "FORM") && df$FORM != -99 && df$FORM == 0) FRELFORM <- 1 + thetas[12]
    FRELSEX <- 1
    if (any(names(df) == "SEX") && df$SEX != -99 && df$SEX == 2) FRELSEX <- 1 + thetas[14]
    FRELGENO4 <- 1
    if (any(names(df) == "GENO4") && df$GENO4 != -99 && df$GENO4 == 1) FRELGENO4 <- 1 + thetas[13]
    FREL <- thetas[1] * FRELSEX * FRELFORM * FRELGENO4
    if (any(names(df) == "WT") && df$WT != -99) {
      TVCL <- thetas[4] * (df$WT / 75)^thetas[2]
    } else {
      TVCL <- thetas[4]
    }
    if (any(names(df) == "GENO1") && df$GENO1 != -99 && df$GENO1 == 1) TVCL <- TVCL * (1 + thetas[8])
    if (any(names(df) == "GENO3") && df$GENO3 != -99 && df$GENO3 == 1) TVCL <- TVCL * (1 + thetas[9])
    if (any(names(df) == "GENO4") && df$GENO4 != -99 && df$GENO4 == 1) TVCL <- TVCL * (1 + thetas[10])
    CL <- CLFOOD * TVCL
    if (any(names(df) == "WT") && df$WT != -99) {
      V <- thetas[5] * (df$WT / 75)^thetas[3]
    } else {
      V <- thetas[5]
    }
    list(CL = CL, FREL = FREL, V = V)
  }

  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("CL", "FREL", "V"), quiet = TRUE)
  )
  generated <- eval(parse(text = out$code))

  dfData <- read.csv(
    system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
  )
  covs <- c("WT", "SEX", "FOOD", "FORM", "GENO1", "GENO3", "GENO4")
  dfCovs <- setupDfCovs(dfData, covariates = covs, idVar = "ID")
  dfCovs <- dfCovs[, setdiff(names(dfCovs), "COVARIATEGROUPS"), drop = FALSE]

  # getForestDFSCM() evaluates the function on an all-missing row when no
  # dfRefRow is given, so that row has to agree too.
  refRow <- dfCovs[1, , drop = FALSE]
  refRow[, ] <- -99
  rows <- rbind(dfCovs, refRow)

  thetas <- run7Thetas()
  for (i in seq_len(nrow(rows))) {
    a <- unlist(handWritten(thetas, rows[i, , drop = FALSE]))
    b <- unlist(generated(thetas, rows[i, , drop = FALSE]))
    expect_equal(b[names(a)], a, tolerance = 1e-12)
  }
})

test_that("the generated function drives getForestDFSCM()", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE)
  )
  fun <- eval(parse(text = out$code))

  dfData <- read.csv(
    system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
  )
  dfCovs <- setupDfCovs(dfData, covariates = c("WT", "FOOD"), idVar = "ID")
  dfSamples <- getSamples(
    system.file("extdata", "SimVal/run7.cov", package = "PMXForest"),
    system.file("extdata", "SimVal/run7.ext", package = "PMXForest"),
    n = 10
  )

  res <- getForestDFSCM(
    dfCovs,
    functionList = list(fun),
    functionListName = out$functionListName,
    noBaseThetas = out$noBaseThetas, dfParameters = dfSamples
  )
  expect_s3_class(res, "data.frame")
  expect_setequal(as.character(unique(res$PARAMETER)), c("CL", "V"))
  expect_true(all(is.finite(res$POINT)))
})

test_that("parameters is validated against what $PK assigns", {
  expect_error(
    suppressWarnings(createParamFunction(modFile, parameters = "NOPE", quiet = TRUE)),
    "Not assigned in the \\$PK block"
  )
})

test_that("the THETA count comes from the .ext file when one is given", {
  out <- suppressWarnings(createParamFunction(
    modFile,
    parameters = "CL", quiet = TRUE,
    extFile = system.file("extdata", "SimVal/run7.ext", package = "PMXForest")
  ))
  expect_equal(out$noBaseThetas, 14)
})

test_that("a THETA index beyond the declared count is refused", {
  f <- tempMod(c("$INPUT ID DV", "$PK", "CL = THETA(9)", "$THETA 1 2"))
  expect_error(
    createParamFunction(f, quiet = TRUE),
    "declares 2 THETA\\(s\\) but \\$PK references THETA\\(9\\)"
  )
})

test_that("a model without $PK is refused", {
  f <- tempMod(c("$INPUT ID DV", "$PRED", "Y = THETA(1) + ETA(1) + EPS(1)"))
  expect_error(createParamFunction(f, quiet = TRUE), "No \\$PK record")
})

test_that("$ERROR is never read, so a model needing an ODE still converts", {
  f <- tempMod(c(
    "$INPUT ID DV WT", "$PK", "CL = THETA(1)", "V = THETA(2)",
    "$ERROR", "CP = A(2)*1000/V", "Y = CP + EPS(1)", "$THETA 1 2"
  ))
  out <- createParamFunction(f, quiet = TRUE)
  expect_equal(out$functionListName, c("CL", "V"))
})

test_that("unsupported syntax inside $PK is refused with a file and line", {
  f <- tempMod(c(
    "$INPUT ID DV", "$PK", "CL = THETA(1)", "CALL MYSUB(CL)",
    "$THETA 1"
  ))
  expect_error(createParamFunction(f, quiet = TRUE), ":4")
})

test_that("the emitted source carries provenance and an extension point", {
  out <- suppressWarnings(createParamFunction(modFile, quiet = TRUE))
  code <- paste(out$code, collapse = "\n")
  expect_match(code, "Generated by PMXForest::createParamFunction")
  expect_match(code, "run7\\.mod:18") # SEX provenance
  expect_match(code, "Secondary parameters: add yours below")
  expect_match(code, "ETA\\(\\) -> 0")
  # one-line IFs stay on one line, so the source diffs against $PK
  expect_true(any(grepl(
    "^\\s*if \\(SEX == 2\\) FRELSEX <- 1 \\+ thetas\\[14\\]$",
    out$code
  )))
})

test_that("file = writes the source and print() renders it", {
  f <- withr::local_tempfile(fileext = ".R")
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", file = f, quiet = TRUE)
  )
  expect_true(file.exists(f))
  expect_equal(readLines(f), as.character(out$code))
  expect_output(print(out$code), "paramFunction <- function")
})

test_that("functionName controls the name of the generated function", {
  out <- suppressWarnings(
    createParamFunction(modFile,
      parameters = "CL", functionName = "myPF",
      quiet = TRUE
    )
  )
  expect_true(any(grepl("^myPF <- function", out$code)))
})

test_that("quiet = FALSE reports the covariates and their references", {
  ## CL reaches five of run7's seven covariates; SEX and FORM belong to FREL
  ## and are pruned away with it.
  expect_message(
    suppressWarnings(createParamFunction(modFile, parameters = "CL")),
    "5 covariate\\(s\\)"
  )
  expect_message(
    suppressWarnings(createParamFunction(modFile)),
    "7 covariate\\(s\\)"
  )
})

test_that("the exponential-IIV map covers the returned parameters only", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE)
  )
  expect_equal(out$etaMap[["CL"]], 3L)
  expect_equal(out$etaMap[["V"]], 4L)
  expect_false("FREL" %in% names(out$etaMap))
})

test_that("a bundled FREM model is refused, however well its $PK would convert", {
  # run22-3's structural block converts perfectly well - that is the trap. Its
  # covariate effects are in OMEGA, so the function it would produce describes
  # none of the covariates the model was built for.
  fremFile <- system.file("extdata", "SimVal/run22-3.mod", package = "PMXForest")
  expect_error(
    createParamFunction(fremFile, parameters = c("CL", "V"), quiet = TRUE),
    "is a FREM model"
  )
})

test_that("MU-referenced structural code converts", {
  # The idiom run22-3 used to cover, without the FREM data item: a parameter
  # defined through the MU layer rather than directly, as SAEM and IMP models
  # normally are.
  mod <- c(
    "$PROBLEM mu", "$INPUT ID TIME DV AMT WT FOOD",
    "$DATA d.csv IGNORE=@", "$PK",
    "IF(FOOD.EQ.1) FOODCL = 1  ; Most common",
    "IF(FOOD.EQ.0) FOODCL = (1 + THETA(3))",
    "TVCL = THETA(1)*(WT/75)**0.75",
    "TVV  = THETA(2)",
    "MU_4 = LOG(TVCL)",
    "MU_5 = LOG(TVV)",
    "CL   = FOODCL * EXP(MU_4 + ETA(4))",
    "V    = EXP(MU_5 + ETA(5))",
    "$THETA (0,7) (0,3) (-1,0.2)"
  )
  out <- createParamFunction(tempMod(mod), parameters = c("CL", "V"), quiet = TRUE)
  expect_equal(out$noBaseThetas, 3)
  expect_true(all(c("WT", "FOOD") %in% names(out$covRef)))
  fun <- eval(parse(text = out$code))
  v <- fun(thetas = c(7, 3, 0.2), df = data.frame(WT = 80, FOOD = 0))
  expect_true(all(is.finite(unlist(v))))
  # MU_4 = log(TVCL), so CL at the reference is FOODCL * TVCL
  ref <- fun(thetas = c(7, 3, 0.2), df = data.frame(WT = 75, FOOD = 1))
  expect_equal(ref$CL, 7)
  expect_equal(ref$V, 3)
})

test_that("missVal is honoured throughout", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", missVal = -999, quiet = TRUE)
  )
  expect_true(any(grepl("!= -999", out$code)))
  fun <- eval(parse(text = out$code))
  thetas <- run7Thetas()
  # -999 marks the covariate inactive, so the reference weight is used
  expect_equal(
    fun(thetas, data.frame(WT = -999, FOOD = -999)),
    fun(thetas, data.frame(WT = 75, FOOD = 1))
  )
})

# ---------------------------------------------------------------------------
# secondary parameters
# ---------------------------------------------------------------------------

test_that("a secondary snippet is spliced in and returned", {
  out <- suppressWarnings(createParamFunction(
    modFile,
    parameters = c("CL", "V"), quiet = TRUE,
    secondary = list(AUC = "df$DOSE / CL", KEL = "CL / V")
  ))
  expect_equal(out$functionListName, c("CL", "V", "AUC", "KEL"))
  expect_equal(out$primaryNames, c("CL", "V"))
  expect_equal(out$secondaryNames, c("AUC", "KEL"))

  code <- paste(out$code, collapse = "\n")
  expect_match(code, "AUC <- local\\(\\{ df\\$DOSE / CL \\}\\)")
  expect_match(code, "KEL <- local\\(\\{ CL / V \\}\\)")
  expect_no_match(code, "add yours below") # extension point replaced

  fun <- eval(parse(text = out$code))
  expect_equal(names(formals(fun)), c("thetas", "df", "..."))
  thetas <- run7Thetas()
  v <- fun(thetas, data.frame(WT = 90, FOOD = 0, DOSE = 160))
  expect_named(v, c("CL", "V", "AUC", "KEL"))
  expect_equal(v$AUC, 160 / v$CL)
  expect_equal(v$KEL, v$CL / v$V)
})

test_that("secondaries are evaluated in order and can use earlier ones", {
  out <- suppressWarnings(createParamFunction(
    modFile,
    parameters = "CL", quiet = TRUE,
    secondary = list(KEL = "CL / 100", HALFLIFE = "log(2) / KEL")
  ))
  fun <- eval(parse(text = out$code))
  v <- fun(run7Thetas(), data.frame(WT = 75, FOOD = 1))
  expect_equal(v$HALFLIFE, log(2) / v$KEL)
})

test_that("a config-list secondary binds its constants ahead of the source", {
  out <- suppressWarnings(createParamFunction(
    modFile,
    parameters = "CL", quiet = TRUE,
    secondary = list(AUC = list(source = "dose / CL", dose = 240))
  ))
  code <- paste(out$code, collapse = "\n")
  expect_match(code, "AUC <- local\\(\\{")
  expect_match(code, "dose <- 240")
  expect_equal(out$functionListName, c("CL", "AUC"))

  fun <- eval(parse(text = out$code))
  v <- fun(run7Thetas(), data.frame(WT = 75, FOOD = 1))
  expect_equal(v$AUC, 240 / v$CL)
})

test_that("a config-list constant can be referenced by a file source", {
  rf <- withr::local_tempfile(fileext = ".R")
  writeLines(c("# uses the injected `tau`", "auc <- dose / CL", "auc / tau"), rf)
  out <- suppressWarnings(createParamFunction(
    modFile,
    parameters = "CL", quiet = TRUE,
    secondary = list(CAVG = list(source = rf, dose = 100, tau = 24))
  ))
  file.remove(rf)
  fun <- eval(parse(text = out$code))
  v <- fun(run7Thetas(), data.frame(WT = 75, FOOD = 1))
  expect_equal(v$CAVG, (100 / v$CL) / 24)
})

test_that("a bare-string secondary still works unchanged (back-compat)", {
  b <- suppressWarnings(createParamFunction(
    modFile,
    parameters = "CL", quiet = TRUE,
    secondary = list(AUC = "80 / CL")
  ))
  l <- suppressWarnings(createParamFunction(
    modFile,
    parameters = "CL", quiet = TRUE,
    secondary = list(AUC = list(source = "80 / CL"))
  )) # source only, no consts
  expect_identical(
    grep("AUC <- local", b$code, value = TRUE),
    grep("AUC <- local", l$code, value = TRUE)
  )
})

test_that("a secondary read from a file is inlined verbatim", {
  rf <- withr::local_tempfile(fileext = ".R")
  writeLines(c(
    "# a small derived quantity",
    "scale <- 1000",
    "scale * V / CL"
  ), rf)
  out <- suppressWarnings(createParamFunction(
    modFile,
    parameters = c("CL", "V"), quiet = TRUE,
    secondary = list(MRT = rf)
  ))
  code <- paste(out$code, collapse = "\n")
  expect_match(code, "MRT <- local\\(\\{")
  expect_match(code, "inlined from ")
  expect_match(code, "a small derived quantity") # the file's own comment
  # the artifact does not depend on the file still being on disk
  file.remove(rf)
  fun <- eval(parse(text = out$code))
  v <- fun(run7Thetas(), data.frame(WT = 75, FOOD = 1))
  expect_equal(v$MRT, 1000 * v$V / v$CL)
})

test_that("the generated function with secondaries drives getForestDFSCM()", {
  out <- suppressWarnings(createParamFunction(
    modFile,
    parameters = c("CL", "V"), quiet = TRUE,
    secondary = list(AUC = "80 / CL")
  ))
  fun <- eval(parse(text = out$code))
  dfData <- read.csv(
    system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
  )
  dfCovs <- setupDfCovs(dfData, covariates = c("WT", "FOOD"), idVar = "ID")
  dfSamples <- getSamples(
    system.file("extdata", "SimVal/run7.cov", package = "PMXForest"),
    system.file("extdata", "SimVal/run7.ext", package = "PMXForest"),
    n = 10
  )
  res <- getForestDFSCM(
    dfCovs,
    functionList = list(fun), functionListName = out$functionListName,
    noBaseThetas = out$noBaseThetas, dfParameters = dfSamples
  )
  expect_setequal(as.character(unique(res$PARAMETER)), c("CL", "V", "AUC"))
  expect_true(all(is.finite(res$POINT)))
})

test_that("quiet = FALSE announces each secondary", {
  expect_message(
    suppressWarnings(createParamFunction(
      modFile,
      parameters = "CL",
      secondary = list(AUC = "80 / CL")
    )),
    "secondary AUC: inline snippet"
  )
})

test_that("a NONMEM-supplied symbol can be pinned through covRef", {
  ## MIXNUM is the subpopulation index of a $MIX model. NONMEM supplies it, so
  ## it is never in $INPUT, and $PK never assigns it - which used to make every
  ## mixture model unusable. It is a legitimate thing to pin, though: fixing it
  ## says "show me subpopulation 1", and varying it across dfCovs rows plots
  ## each subpopulation in turn.
  mod <- c(
    "$PROBLEM mixture", "$INPUT ID TIME DV AMT WT",
    "$DATA data.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)",
    "IF(MIXNUM.EQ.2) TVCL = THETA(2)",
    "CL = TVCL*(WT/75)**0.75",
    "$THETA (0,7) (0,3)", "$MIX"
  )

  ## Refused when nothing says what MIXNUM should be ...
  expect_error(
    createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE),
    "MIXNUM"
  )

  ## ... and accepted when the caller pins it.
  out <- createParamFunction(
    tempMod(mod),
    parameters = "CL", covRef = list(MIXNUM = 1), quiet = TRUE
  )
  expect_true("MIXNUM" %in% names(out$covRef))
  expect_equal(out$covRef$MIXNUM$value, 1)
  expect_equal(out$covRef$MIXNUM$source, "supplied through covRef")

  ## The pinned value reaches the generated code as an ordinary covariate, so
  ## the subpopulation can be chosen per row.
  f <- eval(parse(text = paste(out$code, collapse = "\n")))
  sub1 <- f(thetas = c(7, 3), df = data.frame(WT = 75, MIXNUM = 1))
  sub2 <- f(thetas = c(7, 3), df = data.frame(WT = 75, MIXNUM = 2))
  expect_equal(sub1$CL, 7)
  expect_equal(sub2$CL, 3)
})

test_that("every unbound symbol is reported at once, not one per run", {
  ## A $MIX model typically reads both MIXNUM and MIXEST. Reporting only the
  ## first means the caller pins it, re-runs, is told about the second, pins
  ## that, re-runs... One error should name all of them.
  mod <- c(
    "$PROBLEM mixture", "$INPUT ID TIME DV AMT WT",
    "$DATA data.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)",
    "IF(MIXNUM.EQ.2) TVCL = THETA(2)",
    "IF(MIXEST.EQ.2) TVCL = TVCL*1.1",
    "CL = TVCL*(WT/75)**0.75",
    "$THETA (0,7) (0,3)", "$MIX"
  )
  err <- tryCatch(
    createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "MIXNUM")
  expect_match(err, "MIXEST")
  ## and it hands back something that can be pasted straight in
  expect_match(err, "covRef = list(MIXNUM = <value>, MIXEST = <value>)", fixed = TRUE)

  ## pinning both works
  out <- createParamFunction(
    tempMod(mod),
    parameters = "CL",
    covRef = list(MIXNUM = 1, MIXEST = 1), quiet = TRUE
  )
  expect_true(all(c("MIXNUM", "MIXEST") %in% names(out$covRef)))
})

test_that("symbols to pin and covariates with no reference are reported together", {
  ## Two different routes to the same remedy: MIXNUM/MIXEST are not in $INPUT
  ## at all, while TIME is but has no rule that yields a reference. Both are
  ## fixed by covRef, so both belong in one error with one covRef to paste.
  mod <- c(
    "$PROBLEM mixture", "$INPUT ID TIME DV AMT WT",
    "$DATA data.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)",
    "IF(MIXNUM.EQ.2) TVCL = THETA(2)",
    "IF(MIXEST.EQ.2) TVCL = TVCL*1.1",
    "CL = TVCL*(WT/75)**0.75 + TIME*THETA(3)",
    "$THETA (0,7) (0,3) (0,0.1)", "$MIX"
  )
  err <- tryCatch(
    createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE),
    error = function(e) conditionMessage(e)
  )
  for (nm in c("MIXNUM", "MIXEST", "TIME")) expect_match(err, nm)
  expect_match(
    err,
    "covRef = list(MIXNUM = <value>, MIXEST = <value>, TIME = <value>)",
    fixed = TRUE
  )
  ## and the message says which of them NONMEM supplies, since that is the
  ## part a reader cannot work out from the control stream alone
  expect_match(err, "supplied by NONMEM")

  out <- createParamFunction(
    tempMod(mod),
    parameters = "CL",
    covRef = list(MIXNUM = 1, MIXEST = 1, TIME = 0), quiet = TRUE
  )
  expect_true(all(c("MIXNUM", "MIXEST", "TIME") %in% names(out$covRef)))
})

test_that("a FREM model is refused", {
  ## PsN builds FREM models with a FREMTYPE data item and a marked block of
  ## generated code. Translating their $PK succeeds and is faithful, but it is
  ## not the covariate model the caller is after, so say so rather than letting
  ## them find out from a 57-entry return list.
  frem <- c(
    "$PROBLEM FREM", "$INPUT ID TIME DV AMT WT FREMTYPE",
    "$DATA frem.dta IGNORE=@", "$PK",
    "TVCL = THETA(1)", "CL = TVCL*(WT/75)**0.75",
    "$THETA (0,7)"
  )
  expect_error(
    createParamFunction(tempMod(frem), parameters = "CL", quiet = TRUE),
    "is a FREM model"
  )
  expect_error(
    createParamFunction(tempMod(frem), parameters = "CL", quiet = TRUE),
    "createFREMParamFunction"
  )

  ## and an ordinary model is untouched
  plain <- c(
    "$PROBLEM plain", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)", "CL = TVCL*(WT/75)**0.75",
    "$THETA (0,7)"
  )
  expect_no_error(
    createParamFunction(tempMod(plain), parameters = "CL", quiet = TRUE)
  )
})

test_that("a branch that departs from an unconditional default is not the reference", {
  ## From a real model:
  ##   IND = 0 ; COV=3 or missing (-99)
  ##   IF(COV.EQ.2) IND = 1
  ## IND is an indicator, so its identity is 0 and the unconditional
  ## assignment is the reference state; the IF is the departure from it. Rule 2b
  ## reads "a branch assigning 1" as the reference category, which here names
  ## the treated level - and said so with confidence, silently pinning the
  ## reference subject to the wrong group. There is no way to recover the right
  ## covariate value from $PK (the comment says 3, or missing), so the honest
  ## outcome is a proposal with a warning, not a confident wrong answer.
  mod <- c(
    "$PROBLEM indicator", "$INPUT ID TIME DV AMT COV",
    "$DATA d.csv IGNORE=@", "$PK",
    "IND = 0",
    "IF(COV.EQ.2) IND = 1",
    "TVCL = THETA(1) + IND*THETA(2)",
    "CL = TVCL",
    "$THETA (0,7) (-1,0.1)"
  )
  out <- suppressWarnings(
    createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE)
  )
  expect_false(out$covRef$COV$value == 2)
  expect_false(out$covRef$COV$confident)
  expect_warning(
    createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE),
    "inferred rather than read"
  )

  ## The ordinary SCM shape still resolves confidently: no unconditional
  ## default, so the branch really is the reference category.
  scm <- c(
    "$PROBLEM scm", "$INPUT ID TIME DV AMT FOOD",
    "$DATA d.csv IGNORE=@", "$PK",
    "IF(FOOD.EQ.1) CLFOOD = 1",
    "IF(FOOD.EQ.0) CLFOOD = (1 + THETA(2))",
    "CL = THETA(1)*CLFOOD",
    "$THETA (0,7) (-1,0.1)"
  )
  o2 <- createParamFunction(tempMod(scm), parameters = "CL", quiet = TRUE)
  expect_equal(o2$covRef$FOOD$value, 1)
  expect_true(o2$covRef$FOOD$confident)
})

test_that("verbatim FORTRAN can be ignored on request", {
  ## Verbatim code is refused by default because it can define variables the
  ## rest of the block reads. Often it does not - a solver directive such as
  ## MXSTP01 has nothing to do with the parameter algebra - so the caller can
  ## say so. It stays opt-in: the parser cannot tell the two apart.
  mod <- c(
    "$PROBLEM verbatim", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PK",
    '  "FIRST',
    '  " USE PRDATA, ONLY: MXSTP01',
    '  " MXSTP01=2147483647',
    "TVCL = THETA(1)", "CL = TVCL*(WT/75)**0.75",
    "$THETA (0,7)"
  )
  expect_error(
    createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE),
    "verbatim FORTRAN"
  )
  out <- createParamFunction(
    tempMod(mod),
    parameters = "CL", ignoreVerbatim = TRUE, quiet = TRUE
  )
  f <- eval(parse(text = paste(out$code, collapse = "\n")))
  expect_equal(f(thetas = 7, df = data.frame(WT = 75))$CL, 7)
})

test_that("OMEGA(i,j) resolves from the .ext file", {
  ## A model that builds a correlation by hand needs the OMEGA elements. They
  ## are fixed at the final estimates, so they fold to literals.
  extFile <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")
  mod <- c(
    "$PROBLEM omega", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PK",
    "SD = SQRT(OMEGA(2,2))",
    "CL = THETA(1)*(1 + SD)",
    "$THETA (0,7)"
  )
  out <- createParamFunction(
    tempMod(mod),
    parameters = "CL", extFile = extFile, quiet = TRUE
  )
  f <- eval(parse(text = paste(out$code, collapse = "\n")))
  thetas <- rep(1, out$noBaseThetas)
  expect_equal(f(thetas = thetas, df = data.frame(WT = 75))$CL,
    1 * (1 + sqrt(0.255608)),
    tolerance = 1e-6
  )
  ## symmetric: OMEGA(2,3) is the same element as OMEGA(3,2). run7's value
  ## there is 0, which a broken lookup could also produce, so the value is
  ## pinned against a purpose-built .ext below rather than here.
  mod2 <- sub("OMEGA(2,2)", "OMEGA(2,3)", mod, fixed = TRUE)
  expect_no_error(
    createParamFunction(tempMod(mod2),
      parameters = "CL",
      extFile = extFile, quiet = TRUE
    )
  )
  ## and without an .ext there is nothing to resolve it from
  expect_error(
    createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE),
    "OMEGA"
  )
})

test_that("SIGMA(), and either index order of an off-diagonal, fold correctly", {
  ## The .ext holds the lower triangle only, so OMEGA(1,2) can only be found
  ## by trying OMEGA(2,1). An .ext built here rather than the bundled one:
  ## run7's off-diagonal is 0, so a lookup that silently resolved to the wrong
  ## element - or to nothing - would give the same answer as a correct one.
  ext <- withr::local_tempfile(fileext = ".ext")
  writeLines(c(
    "TABLE NO.     1: First Order Conditional Estimation",
    " ITERATION    THETA1       SIGMA(1,1)   OMEGA(1,1)   OMEGA(2,1)   OMEGA(2,2)",
    "  0.0000E+00  7.0000E+00   1.0000E-01   2.0000E-01   3.0000E-01   4.0000E-01",
    " -1.0000E+09  7.0000E+00   1.1000E-01   2.2000E-01   3.3000E-01   4.4000E-01"
  ), ext)

  build <- function(expr) {
    m <- c(
      "$PROBLEM matrices", "$INPUT ID TIME DV AMT",
      "$DATA d.csv IGNORE=@", "$PK",
      paste0("CL = THETA(1) + ", expr),
      "$THETA (0,7)"
    )
    out <- createParamFunction(tempMod(m),
      parameters = "CL", extFile = ext, quiet = TRUE
    )
    f <- eval(parse(text = paste(out$code, collapse = "\n")))
    f(thetas = rep(0, out$noBaseThetas), df = data.frame(ID = 1))$CL
  }

  ## the final-estimate row, not the initial one
  expect_equal(build("SIGMA(1,1)"), 0.11)
  expect_equal(build("OMEGA(2,2)"), 0.44)
  ## present in the file as OMEGA(2,1); both orders must give that element,
  ## and it is distinct from every other value in the file
  expect_equal(build("OMEGA(2,1)"), 0.33)
  expect_equal(build("OMEGA(1,2)"), 0.33)

  ## an element the .ext does not carry is named, not silently defaulted
  expect_error(build("OMEGA(9,9)"), "OMEGA\\(9,9\\)")
  expect_error(build("SIGMA(1,1) + OMEGA(9,9)"), "\\.ext")
})

test_that("a negated most-common branch says so rather than blaming the IF()s", {
  ## IF(INH.NE.1) INHCOV = 1  ; Most common
  ## The marker is there, on a negation, so the rule that reads a level from an
  ## .EQ. test cannot fire. Reporting that as "level not tested by any IF()"
  ## sends the reader looking for a missing branch that is not the problem.
  mod <- c(
    "$PROBLEM negated", "$INPUT ID TIME DV AMT INH",
    "$DATA d.csv IGNORE=@", "$PK",
    "IF(INH.NE.1) INHCOV = 1  ; Most common",
    "IF(INH.EQ.1) INHCOV = (1 + THETA(2))",
    "CL = THETA(1)*INHCOV",
    "$THETA (0,7) (-1,0.2)"
  )
  out <- suppressWarnings(
    createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE)
  )
  expect_false(out$covRef$INH$confident)
  expect_match(out$covRef$INH$source, "Most common")
  expect_match(out$covRef$INH$source, "negat")
})

test_that("MU-referenced parameters get an etaMap entry", {
  ## EXP(MU_6 + ETA(6)) is EXP(MU_6) * EXP(ETA(6)), so dividing a tabled
  ## individual value by exp(ETA(6)) recovers the typical value exactly as for
  ## the multiplicative idiom. Without this, verifyParamFunction() has nothing
  ## to compare against on any IMP or SAEM model - which is most of them.
  mod <- c(
    "$PROBLEM mu", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)", "TVV = THETA(2)",
    "MU_6 = LOG(TVCL)", "MU_7 = LOG(TVV)",
    "COVEFF = 1",
    "CL = EXP(MU_6 + ETA(6))",
    "V  = COVEFF * EXP(MU_7 + ETA(7))",
    "$THETA (0,7) (0,3)"
  )
  out <- suppressWarnings(
    createParamFunction(tempMod(mod), parameters = c("CL", "V"), quiet = TRUE)
  )
  expect_equal(unname(out$etaMap[["CL"]]), 6L)
  expect_equal(unname(out$etaMap[["V"]]), 7L)
})

test_that("an eta that is scaled inside EXP() is not treated as separable", {
  ## From a real model's hand-built correlation:
  ##   FREL = TVFREL * EXP(CORR*SD_FREL/SD_KA*ETA(3) + SQRT(1-CORR**2)*ETA(2))
  ## Neither eta is a bare additive term, so FREL / exp(ETA(n)) does not give
  ## the typical value and no entry may be claimed.
  mod <- c(
    "$PROBLEM chol", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVFREL = THETA(1)", "CORR = THETA(2)",
    "FREL = TVFREL * EXP(CORR*ETA(3) + ETA(2))",
    "CL = FREL*THETA(3)",
    "$THETA (0,1) (0,0.6) (0,7)"
  )
  out <- suppressWarnings(
    createParamFunction(tempMod(mod), parameters = c("FREL", "CL"), quiet = TRUE)
  )
  expect_false("FREL" %in% names(out$etaMap))

  ## A *single* eta, scaled. The case above has two etas, so the eta-count
  ## guard fires first and the bare-additive-term guard is never reached -
  ## which is the one this test is named for. Without it,
  ## verifyParamFunction() would divide by exp(eta) where exp(0.5*eta) was
  ## meant, and mis-verify silently.
  one <- c(
    "$PROBLEM scaled", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)", "SC = THETA(2)",
    "CL = TVCL * EXP(SC*ETA(3))",
    "$THETA (0,7) (0,0.5)"
  )
  o2 <- suppressWarnings(
    createParamFunction(tempMod(one), parameters = "CL", quiet = TRUE)
  )
  expect_false("CL" %in% names(o2$etaMap))
})

## --- pruning to the requested parameters ---------------------------------

test_that("only the statements a requested parameter needs are emitted", {
  ## run7's V depends on TVV <- VCOV1 <- VWT <- WT, and on nothing else. The
  ## FREL and CL machinery, and MAT/D1/KA/S2, have no bearing on it.
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "V", quiet = TRUE)
  )
  code <- paste(out$code, collapse = "\n")

  for (kept in c("VWT", "VCOV1", "TVV")) {
    expect_match(code, kept, info = kept)
  }
  for (dropped in c("FRELCOV", "FRELSEX", "CLFOOD", "CLGENO1", "TVMAT", "TVD1", "KA")) {
    expect_false(grepl(dropped, code, fixed = TRUE), info = dropped)
  }
})

test_that("pruning drops the covariates the parameter does not use", {
  ## V reaches only WT. CL reaches WT, FOOD and the three genotype dummies,
  ## but neither SEX nor FORM, which belong to FREL.
  v <- suppressWarnings(createParamFunction(modFile, parameters = "V", quiet = TRUE))
  expect_equal(names(v$covRef), "WT")

  cl <- suppressWarnings(createParamFunction(modFile, parameters = "CL", quiet = TRUE))
  expect_setequal(names(cl$covRef), c("WT", "FOOD", "GENO1", "GENO3", "GENO4"))
  expect_false(any(c("SEX", "FORM") %in% names(cl$covRef)))
})

test_that("a pruned function returns the same numbers as an unpruned one", {
  ## Pruning may only remove statements that cannot affect the result.
  full <- suppressWarnings(createParamFunction(modFile, quiet = TRUE))
  part <- suppressWarnings(createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE))
  thetas <- run7Thetas()
  fFull <- eval(parse(text = paste(full$code, collapse = "\n")))
  fPart <- eval(parse(text = paste(part$code, collapse = "\n")))

  row <- data.frame(
    WT = 84, SEX = 2, FOOD = 0, FORM = 0,
    GENO1 = 1, GENO3 = 0, GENO4 = 0
  )
  a <- fFull(thetas = thetas, df = row)
  b <- fPart(thetas = thetas, df = row)
  expect_equal(b$CL, a$CL)
  expect_equal(b$V, a$V)
  expect_setequal(names(b), c("CL", "V"))
})

test_that("a re-assignment chain survives pruning", {
  ## TVCL is assigned twice - THETA(4)*CLCOV1, then CLCOV*TVCL. Keeping only
  ## the last would silently drop the covariate effect.
  out <- suppressWarnings(createParamFunction(modFile, parameters = "CL", quiet = TRUE))
  code <- paste(out$code, collapse = "\n")
  expect_match(code, "CLCOV1", fixed = TRUE)
  expect_match(code, "CLCOV <-", fixed = TRUE)
  expect_equal(sum(grepl("^\\s*TVCL <-", out$code)), 2L)
})

test_that("an IF block is kept whole, and its condition becomes a dependency", {
  mod <- c(
    "$PROBLEM prune", "$INPUT ID TIME DV AMT WT SEX",
    "$DATA d.csv IGNORE=@", "$PK",
    "IF(SEX.EQ.1) SEXCL = 1  ; Most common",
    "IF(SEX.EQ.2) SEXCL = (1 + THETA(2))",
    "OTHER = THETA(3)*WT",
    "CL = THETA(1)*SEXCL",
    "V  = OTHER",
    "$THETA (0,7) (-1,0.2) (0,1)"
  )
  out <- createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE)
  code <- paste(out$code, collapse = "\n")
  expect_match(code, "SEXCL", fixed = TRUE)
  expect_false(grepl("OTHER", code, fixed = TRUE))
  ## SEX is only ever read in the condition, and is still a covariate
  expect_equal(names(out$covRef), "SEX")

  ## A multi-statement IF/ELSE, so the walk really does descend into `else_`
  ## and not only into a one-line `then`. Losing a branch would drop either an
  ## assignment the parameter needs or a dependency it reads - silently.
  blk <- c(
    "$PROBLEM block", "$INPUT ID TIME DV AMT WT FOOD",
    "$DATA d.csv IGNORE=@", "$PK",
    "IF(FOOD.EQ.1) THEN",
    "  CLF = 1",
    "  SPARE = 99",
    "ELSE",
    "  CLF = (1 + THETA(2)*(WT/75))",
    "  SPARE = 98",
    "ENDIF",
    "CL = THETA(1)*CLF",
    "$THETA (0,7) (-1,0.01)"
  )
  o3 <- suppressWarnings(
    createParamFunction(tempMod(blk), parameters = "CL", quiet = TRUE)
  )
  ## WT is read only inside the ELSE branch, so it is a covariate only if the
  ## walk went there
  expect_true(all(c("FOOD", "WT") %in% names(o3$covRef)))
  expect_match(paste(o3$code, collapse = "\n"), "CLF", fixed = TRUE)
})

test_that("an assignment that cannot reach the parameter is dropped", {
  mod <- c(
    "$PROBLEM dead", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)",
    "CL   = TVCL*(WT/75)**0.75",
    "TVCL = 999",
    "$THETA (0,7)"
  )
  out <- createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE)
  expect_false(any(grepl("999", out$code, fixed = TRUE)))
})

test_that("without `parameters` nothing is pruned", {
  out <- suppressWarnings(createParamFunction(modFile, quiet = TRUE))
  code <- paste(out$code, collapse = "\n")
  for (nm in c("FRELCOV", "CLGENO1", "TVMAT", "TVD1", "S2")) {
    expect_match(code, nm, info = nm)
  }
})

test_that("pruning keeps what a secondary expression reads", {
  ## A secondary block is spliced in after the $PK translation and sees every
  ## structural parameter by name - including intermediates the caller did not
  ## ask for. Pruning to `parameters` alone removed them, and the generated
  ## function then failed at call time with "object 'FREL' not found".
  out <- suppressWarnings(createParamFunction(
    modFile,
    parameters = c("CL", "V"),
    secondary = list(AUC = "80 / (CL / FREL)"),
    quiet = TRUE
  ))
  code <- paste(out$code, collapse = "\n")
  expect_match(code, "FREL", fixed = TRUE)

  row <- data.frame(
    WT = 80, SEX = 1, FOOD = 1, FORM = 1, GENO1 = 0, GENO3 = 0, GENO4 = 0
  )
  v <- eval(parse(text = code))(thetas = run7Thetas(), df = row)
  expect_true(is.finite(v$AUC))

  ## and it is the right number: FREL is recoverable by asking for it, so the
  ## secondary can be checked against the same arithmetic done here
  both <- suppressWarnings(createParamFunction(
    modFile,
    parameters = c("CL", "FREL"), quiet = TRUE
  ))
  w <- eval(parse(text = paste(both$code, collapse = "\n")))(
    thetas = run7Thetas(), df = row
  )
  expect_equal(v$AUC, 80 / (w$CL / w$FREL))

  ## FREL is kept because the secondary needs it, not returned in its own right
  expect_setequal(out$functionListName, c("CL", "V", "AUC"))
})

test_that("a secondary reading nothing extra still prunes", {
  out <- suppressWarnings(createParamFunction(
    modFile,
    parameters = "V", secondary = list(HALF = "V / 2"), quiet = TRUE
  ))
  code <- paste(out$code, collapse = "\n")
  expect_false(grepl("FRELCOV", code, fixed = TRUE))
  expect_setequal(out$functionListName, c("V", "HALF"))
})

test_that("an IF/THEN/ELSE covariate keeps its confident reference", {
  ## The guard that stops a *departure* branch being read as the reference
  ## must not catch the ordinary IF/THEN/ELSE coding, where the ELSE body is
  ## a branch like any other. It did: nmFlatten() records an ELSE body with a
  ## NULL condition, so it looked unconditional, and the rule was skipped.
  ##
  ## The consequence was not a missing reference but a wrong one. SEX would be
  ## proposed as 0, which takes the ELSE branch, so the reference row of the
  ## plot - the denominator of every ratio - is computed for the treated
  ## category.
  mod <- c(
    "$PROBLEM ifelse", "$INPUT ID TIME DV AMT WT SEX",
    "$DATA d.csv IGNORE=@", "$PK",
    "IF(SEX.EQ.1) THEN",
    "  FSEX = 1",
    "ELSE",
    "  FSEX = 1 + THETA(2)",
    "ENDIF",
    "CL = THETA(1)*FSEX*(WT/75)**0.75",
    "$THETA (0,7) (-1,0.3)"
  )
  out <- createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE)
  expect_equal(out$covRef$SEX$value, 1)
  expect_true(out$covRef$SEX$confident)
  expect_match(out$covRef$SEX$source, "identity value")

  ## and the reference really is the untreated category
  f <- eval(parse(text = paste(out$code, collapse = "\n")))
  expect_equal(f(thetas = c(7, 0.3), df = data.frame(WT = 75, SEX = 1))$CL, 7)

  ## the indicator case it was written for still falls through
  ind <- c(
    "$PROBLEM indicator", "$INPUT ID TIME DV AMT COV",
    "$DATA d.csv IGNORE=@", "$PK",
    "IND = 0", "IF(COV.EQ.2) IND = 1",
    "CL  = THETA(1) + IND*THETA(2)",
    "$THETA (0,7) (-1,0.5)"
  )
  o2 <- suppressWarnings(
    createParamFunction(tempMod(ind), parameters = "CL", quiet = TRUE)
  )
  expect_false(o2$covRef$COV$confident)
})

test_that("an eta hidden under a unary minus is still counted", {
  ## etasIn() walks lhs/rhs/cond/args. The parser builds unary minus as
  ## list(type = "unop", op, arg = ...), so an eta under `arg` was invisible,
  ## the one-eta guard passed, and a non-separable expression was claimed.
  mod <- c(
    "$PROBLEM unary", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)", "MU_1 = LOG(TVCL)",
    "CL = EXP(MU_1 + ETA(1) + (-ETA(2)))",
    "$THETA (0,7)"
  )
  out <- suppressWarnings(
    createParamFunction(tempMod(mod), parameters = "CL", quiet = TRUE)
  )
  expect_false("CL" %in% names(out$etaMap))
})

test_that("an eta reaching the exponent through a symbol is not separable", {
  ## The IOV idiom: IOV is assigned from ETA() in an occasion block, then added
  ## inside the exponent. etasIn() counts syntactic ETA() nodes, so IOV
  ## contributed nothing, the one-eta guard passed, and the entry was claimed -
  ## after which CL / exp(ETA1) is TVCL * exp(IOV), not the typical value.
  mod <- c(
    "$PROBLEM iov", "$INPUT ID TIME DV AMT WT OCC",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)",
    "IOV = 0",
    "IF(OCC.EQ.1) IOV = ETA(5)",
    "IF(OCC.EQ.2) IOV = ETA(6)",
    "CL = TVCL*EXP(ETA(1) + IOV)",
    "$THETA (0,7)"
  )
  out <- suppressWarnings(
    createParamFunction(tempMod(mod), parameters = "CL", covRef = list(OCC = 1), quiet = TRUE)
  )
  expect_false("CL" %in% names(out$etaMap))

  ## a plain symbol that carries no randomness is still fine
  ok <- c(
    "$PROBLEM plain", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)", "MU_1 = LOG(TVCL)", "SHIFT = THETA(2)",
    "CL = EXP(MU_1 + SHIFT + ETA(1))",
    "$THETA (0,7) (-1,0.1)"
  )
  o2 <- suppressWarnings(
    createParamFunction(tempMod(ok), parameters = "CL", quiet = TRUE)
  )
  expect_equal(unname(o2$etaMap[["CL"]]), 1L)
})
