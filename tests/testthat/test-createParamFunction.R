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
  expect_named(out, c("code", "functionListName", "primaryNames",
                      "secondaryNames", "noBaseThetas", "covRef",
                      "etaMap", "modFile", "missVal"))
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

test_that("covariates are those in $INPUT that $PK never assigns", {
  out <- suppressWarnings(createParamFunction(modFile, quiet = TRUE))
  expect_setequal(names(out$covRef),
                  c("SEX", "GENO4", "FORM", "FOOD", "WT", "GENO1", "GENO3"))
  # GENO2 is assigned a constant in $PK, so it is a local, not a covariate
  expect_false("GENO2" %in% names(out$covRef))
  # TVCL is assigned in $PK and is never a covariate
  expect_false("TVCL" %in% names(out$covRef))
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
  expect_warning(createParamFunction(modFile, quiet = TRUE),
                 "inferred rather than read")
})

test_that("covRef overrides a derived reference and is recorded as supplied", {
  out <- suppressWarnings(
    createParamFunction(modFile, covRef = list(WT = 70), quiet = TRUE)
  )
  expect_equal(out$covRef$WT$value, 70)
  expect_match(out$covRef$WT$source, "covRef")
  expect_true(any(grepl("else 70", out$code)))
})

test_that("a covariate with no derivable reference is refused, not guessed", {
  # EXPO enters tte_weibull.mod as THETA(4)*EXPO: no branch, no normalisation.
  expect_error(createParamFunction(tteFile, quiet = TRUE),
               "No reference value could be derived")
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
  covs   <- c("WT", "SEX", "FOOD", "FORM", "GENO1", "GENO3", "GENO4")
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
    dfCovs, functionList = list(fun),
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
    modFile, parameters = "CL", quiet = TRUE,
    extFile = system.file("extdata", "SimVal/run7.ext", package = "PMXForest")
  ))
  expect_equal(out$noBaseThetas, 14)
})

test_that("a THETA index beyond the declared count is refused", {
  f <- tempMod(c("$INPUT ID DV", "$PK", "CL = THETA(9)", "$THETA 1 2"))
  expect_error(createParamFunction(f, quiet = TRUE),
               "declares 2 THETA\\(s\\) but \\$PK references THETA\\(9\\)")
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
  f <- tempMod(c("$INPUT ID DV", "$PK", "CL = THETA(1)", "CALL MYSUB(CL)",
                 "$THETA 1"))
  expect_error(createParamFunction(f, quiet = TRUE), ":4")
})

test_that("the emitted source carries provenance and an extension point", {
  out <- suppressWarnings(createParamFunction(modFile, quiet = TRUE))
  code <- paste(out$code, collapse = "\n")
  expect_match(code, "Generated by PMXForest::createParamFunction")
  expect_match(code, "run7\\.mod:18")                    # SEX provenance
  expect_match(code, "Secondary parameters: add yours below")
  expect_match(code, "ETA\\(\\) -> 0")
  # one-line IFs stay on one line, so the source diffs against $PK
  expect_true(any(grepl("^\\s*if \\(SEX == 2\\) FRELSEX <- 1 \\+ thetas\\[14\\]$",
                        out$code)))
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
    createParamFunction(modFile, parameters = "CL", functionName = "myPF",
                        quiet = TRUE)
  )
  expect_true(any(grepl("^myPF <- function", out$code)))
})

test_that("quiet = FALSE reports the covariates and their references", {
  expect_message(
    suppressWarnings(createParamFunction(modFile, parameters = "CL")),
    "7 covariate\\(s\\)"
  )
})

test_that("the exponential-IIV map covers the returned parameters only", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = c("CL", "V"), quiet = TRUE)
  )
  expect_equal(out$etaMap[["CL"]], 3L)
  expect_equal(out$etaMap[["V"]],  4L)
  expect_false("FREL" %in% names(out$etaMap))
})

test_that("the structural $PK of a FREM model converts", {
  # The FREM covariate effects live in OMEGA and are out of scope, but the
  # structural block, including its MU-referencing lines, must not be refused.
  fremFile <- system.file("extdata", "SimVal/run22-3.mod", package = "PMXForest")
  out <- suppressWarnings(
    createParamFunction(fremFile, parameters = c("CL", "V"), quiet = TRUE)
  )
  expect_equal(out$noBaseThetas, 24)
  expect_true(all(c("WT", "FOOD", "FORM") %in% names(out$covRef)))
  fun <- eval(parse(text = out$code))
  v <- fun(thetas = rep(1, 24), df = data.frame(WT = 80, FOOD = 0))
  expect_true(all(is.finite(unlist(v))))
})

test_that("missVal is honoured throughout", {
  out <- suppressWarnings(
    createParamFunction(modFile, parameters = "CL", missVal = -999, quiet = TRUE)
  )
  expect_true(any(grepl("!= -999", out$code)))
  fun <- eval(parse(text = out$code))
  thetas <- run7Thetas()
  # -999 marks the covariate inactive, so the reference weight is used
  expect_equal(fun(thetas, data.frame(WT = -999, FOOD = -999)),
               fun(thetas, data.frame(WT = 75,   FOOD = 1)))
})

# ---------------------------------------------------------------------------
# secondary parameters
# ---------------------------------------------------------------------------

test_that("a secondary snippet is spliced in and returned", {
  out <- suppressWarnings(createParamFunction(
    modFile, parameters = c("CL", "V"), quiet = TRUE,
    secondary = list(AUC = "df$DOSE / CL", KEL = "CL / V")
  ))
  expect_equal(out$functionListName, c("CL", "V", "AUC", "KEL"))
  expect_equal(out$primaryNames,   c("CL", "V"))
  expect_equal(out$secondaryNames, c("AUC", "KEL"))

  code <- paste(out$code, collapse = "\n")
  expect_match(code, "AUC <- local\\(\\{ df\\$DOSE / CL \\}\\)")
  expect_match(code, "KEL <- local\\(\\{ CL / V \\}\\)")
  expect_no_match(code, "add yours below")   # extension point replaced

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
    modFile, parameters = "CL", quiet = TRUE,
    secondary = list(KEL = "CL / 100", HALFLIFE = "log(2) / KEL")
  ))
  fun <- eval(parse(text = out$code))
  v <- fun(run7Thetas(), data.frame(WT = 75, FOOD = 1))
  expect_equal(v$HALFLIFE, log(2) / v$KEL)
})

test_that("a secondary read from a file is inlined verbatim", {
  rf <- withr::local_tempfile(fileext = ".R")
  writeLines(c("# a small derived quantity",
               "scale <- 1000",
               "scale * V / CL"), rf)
  out <- suppressWarnings(createParamFunction(
    modFile, parameters = c("CL", "V"), quiet = TRUE,
    secondary = list(MRT = rf)
  ))
  code <- paste(out$code, collapse = "\n")
  expect_match(code, "MRT <- local\\(\\{")
  expect_match(code, "inlined from ")
  expect_match(code, "a small derived quantity")   # the file's own comment
  # the artifact does not depend on the file still being on disk
  file.remove(rf)
  fun <- eval(parse(text = out$code))
  v <- fun(run7Thetas(), data.frame(WT = 75, FOOD = 1))
  expect_equal(v$MRT, 1000 * v$V / v$CL)
})

test_that("the generated function with secondaries drives getForestDFSCM()", {
  out <- suppressWarnings(createParamFunction(
    modFile, parameters = c("CL", "V"), quiet = TRUE,
    secondary = list(AUC = "80 / CL")
  ))
  fun <- eval(parse(text = out$code))
  dfData <- read.csv(
    system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
  )
  dfCovs    <- setupDfCovs(dfData, covariates = c("WT", "FOOD"), idVar = "ID")
  dfSamples <- getSamples(
    system.file("extdata", "SimVal/run7.cov", package = "PMXForest"),
    system.file("extdata", "SimVal/run7.ext", package = "PMXForest"), n = 10
  )
  res <- getForestDFSCM(
    dfCovs, functionList = list(fun), functionListName = out$functionListName,
    noBaseThetas = out$noBaseThetas, dfParameters = dfSamples
  )
  expect_setequal(as.character(unique(res$PARAMETER)), c("CL", "V", "AUC"))
  expect_true(all(is.finite(res$POINT)))
})

test_that("quiet = FALSE announces each secondary", {
  expect_message(
    suppressWarnings(createParamFunction(
      modFile, parameters = "CL",
      secondary = list(AUC = "80 / CL")
    )),
    "secondary AUC: inline snippet"
  )
})
