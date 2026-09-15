## $PRED models. The block is the same abbreviated-code language as $PK, but it
## computes the prediction and the residual error in the same lines - so `Y` and
## `EPS()` are part of it, and "which assignments are parameters" is not
## answerable from the block alone.

tempMod <- function(lines) {
  f <- withr::local_tempfile(fileext = ".mod", .local_envir = parent.frame())
  writeLines(lines, f)
  f
}

## The same model written both ways. A $PRED block computing CL and V by hand,
## and the ADVAN equivalent whose $PK assigns exactly the same two parameters.
predLines <- c(
  "$PROBLEM one-compartment, hand coded",
  "$INPUT ID TIME DV AMT WT SEX",
  "$DATA d.csv IGNORE=@",
  "$PRED",
  "TVCL = THETA(1)*(WT/70)**THETA(3)",
  "IF(SEX.EQ.1) TVCL = TVCL*(1 + THETA(4))",
  "TVV  = THETA(2)",
  "CL   = TVCL*EXP(ETA(1))",
  "V    = TVV*EXP(ETA(2))",
  "K    = CL/V",
  "F    = AMT/V*EXP(-K*TIME)",
  "Y    = F + F*EPS(1)",
  "$THETA (0,7) (0,30) (0,0.75) (-1,0.2)"
)
pkLines <- sub("^\\$PRED$", "$PK", predLines[!grepl("^(F|Y) ", predLines)])

test_that("a $PRED model translates and gives the same typical values as $PK", {
  p <- createParamFunction(tempMod(predLines),
    parameters = c("CL", "V"), covRef = list(SEX = 0), quiet = TRUE
  )
  k <- createParamFunction(tempMod(pkLines),
    parameters = c("CL", "V"), covRef = list(SEX = 0), quiet = TRUE
  )
  fp <- eval(parse(text = paste(p$code, collapse = "\n")))
  fk <- eval(parse(text = paste(k$code, collapse = "\n")))

  th <- c(7, 30, 0.75, 0.2)
  for (df in list(
    data.frame(WT = 70, SEX = 0), data.frame(WT = 90, SEX = 1),
    data.frame(WT = 55, SEX = 0)
  )) {
    expect_equal(unlist(fp(thetas = th, df = df)), unlist(fk(thetas = th, df = df)))
  }
  ## and the value is the one the algebra says, not merely self-consistent
  expect_equal(fp(thetas = th, df = data.frame(WT = 70, SEX = 0))$CL, 7)
  expect_equal(fp(thetas = th, df = data.frame(WT = 70, SEX = 1))$CL, 7 * 1.2)
})

test_that("$PRED requires `parameters`, and says why", {
  e <- tryCatch(
    createParamFunction(tempMod(predLines), covRef = list(SEX = 0), quiet = TRUE),
    error = conditionMessage
  )
  expect_match(e, "\\$PRED")
  expect_match(e, "parameters")
  ## the reason, not just the requirement - a $PRED block has no structural
  ## signal for which assignments are parameters
  expect_match(e, "residual error|prediction")

  ## $PK is unaffected: NULL still means every assigned variable
  expect_no_error(
    createParamFunction(tempMod(pkLines), covRef = list(SEX = 0), quiet = TRUE)
  )
})

test_that("EPS() and ERR() fold to zero in $PRED, and Y is just an assignment", {
  pins <- list(SEX = 0, AMT = 100, TIME = 0)
  out <- createParamFunction(tempMod(predLines),
    parameters = "Y", covRef = pins, quiet = TRUE
  )
  f <- eval(parse(text = paste(out$code, collapse = "\n")))
  th <- c(7, 30, 0.75, 0.2)
  got <- f(thetas = th, df = data.frame(WT = 70, SEX = 0, AMT = 100, TIME = 0))$Y
  ## Y = F + F*EPS(1) with EPS -> 0 is F = AMT/V at TIME 0
  expect_equal(got, 100 / 30)

  ## ERR() is the same thing under NM-TRAN's other spelling
  errLines <- sub("EPS(1)", "ERR(1)", predLines, fixed = TRUE)
  o2 <- createParamFunction(tempMod(errLines),
    parameters = "Y", covRef = pins, quiet = TRUE
  )
  f2 <- eval(parse(text = paste(o2$code, collapse = "\n")))
  expect_equal(
    f2(thetas = th, df = data.frame(WT = 70, SEX = 0, AMT = 100, TIME = 0))$Y,
    got
  )
})

test_that("EPS() and ERR() are still refused in $PK", {
  for (sym in c("EPS", "ERR")) {
    f <- tempMod(c(
      "$PROBLEM p", "$INPUT ID DV WT", "$DATA d.csv IGNORE=@", "$PK",
      "CL = THETA(1)", paste0("V = THETA(2) + ", sym, "(1)"), "$THETA 1 2"
    ))
    e <- tryCatch(createParamFunction(f, quiet = TRUE), error = conditionMessage)
    expect_match(e, paste0(sym, "\\(\\) cannot appear in \\$PK"))
  }
})

test_that("an eta sharing the exponent with an EPS is not claimed separable", {
  ## etasIn() counts one ETA() here and the bare-eta term is found, so without
  ## an explicit check the entry is claimed - and CL / exp(eta1) still carries
  ## exp(EPS(1)).
  out <- createParamFunction(tempMod(c(
    "$PROBLEM eps-in-exponent", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PRED",
    "TVCL = THETA(1)",
    "CL = TVCL*EXP(ETA(1) + EPS(1))",
    "V  = THETA(2)*EXP(ETA(2))",
    "Y  = CL + EPS(2)",
    "$THETA (0,7) (0,30)"
  )), parameters = c("CL", "V"), quiet = TRUE)

  expect_false("CL" %in% names(out$etaMap))
  ## the ordinary one alongside it still is, so this is not a blanket refusal
  expect_equal(unname(out$etaMap[["V"]]), 2L)
})

test_that("error messages name $PRED, not $PK", {
  ## A refusal that names the wrong record sends the reader to a block that is
  ## not in their model.
  e <- tryCatch(
    createParamFunction(tempMod(c(
      "$PROBLEM p", "$INPUT ID DV", "$DATA d.csv IGNORE=@", "$PRED",
      "CL = THETA(1)", "CALL SOMETHING(X)", "Y = CL", "$THETA 1"
    )), parameters = "CL", quiet = TRUE),
    error = conditionMessage
  )
  expect_match(e, "\\$PRED")
  expect_false(grepl("\\$PK", e))

  ## and the read-before-assign report, which names the block twice
  e2 <- tryCatch(
    createParamFunction(tempMod(c(
      "$PROBLEM p", "$INPUT ID DV WT", "$DATA d.csv IGNORE=@", "$PRED",
      "CL = TVCL*(WT/75)", "TVCL = THETA(1)", "Y = CL", "$THETA 1"
    )), parameters = "CL", quiet = TRUE),
    error = conditionMessage
  )
  expect_match(e2, "\\$PRED")
  expect_false(grepl("\\$PK", e2))
})

test_that("a model with neither $PK nor $PRED names both", {
  e <- tryCatch(
    createParamFunction(tempMod(c(
      "$PROBLEM p", "$INPUT ID DV", "$DATA d.csv IGNORE=@", "$THETA 1"
    )), quiet = TRUE),
    error = conditionMessage
  )
  expect_match(e, "\\$PK")
  expect_match(e, "\\$PRED")
})

test_that("the raw $PRED tree survives nmDeparse(), which PMXFrem relies on", {
  ## nmParsePK() returns the unfolded tree so a downstream emitter can decide
  ## what to do with ETA()/EPS(). An unhandled node type surfaces there as
  ## "Internal error: cannot deparse node type", not as a diagnosis.
  p <- nmParsePK(tempMod(predLines),
    parameters = c("CL", "V", "Y"),
    covRef = list(SEX = 0, AMT = 100, TIME = 0)
  )
  asgn <- Filter(function(s) s$type == "assign", p$statements)
  expect_no_error(src <- vapply(asgn, function(s) nmDeparse(s$rhs), ""))

  ## EPS folds to epsValue, default "0", exactly as ETA folds to etaValue -
  ## so a typical-value emitter needs no special case.
  yRhs <- asgn[[which(vapply(asgn, function(s) s$lhs, "") == "Y")]]$rhs
  expect_equal(nmDeparse(yRhs), "F + F * 0")
  ## and an emitter that wants the residual error back can have it
  expect_equal(nmDeparse(yRhs, epsValue = "EPS1"), "F + F * EPS1")
})

test_that("setupDfRefRow() reads references out of a $PRED model too", {
  ## refModelValues() is a second, independent $PK-only reader, reached by
  ## contRef = "model". Without extending it, this fails on a model
  ## createParamFunction() has just handled.
  f <- tempMod(c(
    "$PROBLEM p", "$INPUT ID TIME DV AMT WT",
    "$DATA d.csv IGNORE=@", "$PRED",
    "CL = THETA(1)*(WT/75)",
    "Y  = CL + EPS(1)",
    "$THETA (0,7)"
  ))
  expect_no_error(
    r <- setupDfRefRow(
      dfCovs = data.frame(WT = -99),
      data = data.frame(ID = 1:20, WT = seq(50, 88, length.out = 20)),
      covariates = "WT", model = f, contRef = "model"
    )
  )
  ## 75 is the normalisation constant the $PRED block divides by, not the
  ## median of the data (69) - so the reference really came out of the model.
  expect_equal(r$WT, 75)
})

test_that("an IOV parameter gets a tvMap entry, so it can still be verified", {
  ## The traditional IOV coding puts three ETA()s in one exponent. None of them
  ## separates, so there is rightly no etaMap entry - but the typical value is
  ## still TVCL, and without a tvMap entry verifyParamFunction() has nothing to
  ## compare a TVCL column against on any IOV model.
  out <- createParamFunction(tempMod(c(
    "$PROBLEM iov", "$INPUT ID TIME DV AMT WT OCC",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVCL = THETA(1)*(WT/70)**THETA(3)",
    "TVV  = THETA(2)",
    "CL = TVCL*EXP(ETA(1) + ETA(2)*(1-OCC) + ETA(3)*OCC)",
    "V  = TVV*EXP(ETA(4))",
    "$THETA (0,7) (0,30) (0,0.75)"
  )), parameters = c("CL", "V"), quiet = TRUE)

  ## OCC is no longer a covariate at all: the whole IOV term folds away, so
  ## the occasion cannot affect the typical value and nothing needs a
  ## reference for it. Before, every IOV model demanded one.
  expect_false("OCC" %in% names(out$covRef))

  expect_equal(unname(out$tvMap[["CL"]]), "TVCL")
  expect_equal(unname(out$tvMap[["V"]]), "TVV")
  ## and CL is still not claimed separable
  expect_false("CL" %in% names(out$etaMap))
  expect_equal(unname(out$etaMap[["V"]]), 4L)

  ## the emitted code collapses rather than carrying 0 * (1 - OCC) + 0 * OCC
  clLine <- grep("^  CL <-", out$code, value = TRUE)
  expect_match(clLine, "CL <- TVCL")
  expect_false(grepl("0 \\* ", clLine))

  ## and the value is unchanged by the collapse
  f <- eval(parse(text = paste(out$code, collapse = "\n")))
  expect_equal(f(thetas = c(7, 30, 0.75), df = data.frame(WT = 70))$CL, 7)
})

test_that("a parameter computed from another still gets no tvMap entry", {
  ## The widening must not turn every assignment into a typical value.
  ## run7 has MAT = TVMAT*EXP(ETA(5)) and then D1 = MAT*(1-TVD1), where TVD1 is
  ## a dimensionless fraction rather than D1's typical value.
  out <- createParamFunction(tempMod(c(
    "$PROBLEM tvmap", "$INPUT ID TIME DV AMT",
    "$DATA d.csv IGNORE=@", "$PK",
    "TVMAT = THETA(1)", "TVD1 = THETA(2)",
    "MAT = TVMAT*EXP(ETA(1))",
    "D1  = MAT*(1 - TVD1)",
    "$THETA (0,2) (0,0.8)"
  )), parameters = c("MAT", "D1"), quiet = TRUE)

  expect_equal(unname(out$tvMap[["MAT"]]), "TVMAT")
  expect_false("D1" %in% names(out$tvMap))
})

test_that("the walkers that ignore an eps node do so deliberately", {
  ## A new node type is skipped silently by every walker whose switch() has a
  ## NULL default. That is the right answer for EPS - it contributes no symbol,
  ## no THETA and no ETA - but it is currently luck rather than intent, and a
  ## walker that started mis-handling it would fail somewhere far away.
  eps <- list(type = "eps", index = 1L)
  expr <- list(
    type = "binop", op = "+",
    lhs = list(type = "sym", name = "F"), rhs = eps
  )
  expect_equal(nmExprSyms(expr), "F")
  expect_equal(PMXForest:::nmMaxTheta(list(
    list(type = "assign", lhs = "Y", rhs = expr)
  )), 0)
  expect_false(PMXForest:::nmHasEta(eps))
  ## and it folds to zero rather than being carried through
  expect_equal(PMXForest:::nmSimplify(eps), list(type = "num", value = 0))
  ## nmDeparse is the one walker that raises on an unknown type, so it needs
  ## the case explicitly - this is the path PMXFrem takes on the raw tree
  expect_equal(nmDeparse(eps), "0")
})
